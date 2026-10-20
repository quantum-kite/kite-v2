"""Run the standalone regression scripts under pytest.

These predate the suite and are written as scripts with a `__main__` block
rather than as `test_` functions, so pytest finds nothing to run in them even
once it collects the files. Five are invoked by name in CI; the two that need a
KITEx run were left out and so ran nowhere. Executing them here puts all seven
behind the single `pytest tests/` entry point without rewriting them.

Each is run as a subprocess, the way CI runs them, so a script that depends on
module-level state or on its own `__main__` behaves identically.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# Self-contained: each builds whatever it needs.
SELF_CONTAINED = [
    "sigma_y_regression.py",
    "operator_ldos_regression.py",
    "operator_ldos_analytic_stochastic_regression.py",
    "operator_dimension_regression.py",
    "hamiltonian_k_convention_regression.py",
]

# Need an example's HDF5 output, which means a prior KITEx run.
NEEDS_KITEX_OUTPUT = {
    "custom_two_mod4_regression.py": "examples/kane_mele_spin_hall-output.h5",
    "custom_two_same_vertex_sign_regression.py":
        "tests/fixtures/graphene_custom_two_xx_check-output.h5",
}


def _run(script):
    env = {"PYTHONPATH": str(REPO / "python"), "MPLBACKEND": "Agg"}
    import os
    env = {**os.environ, **env}
    return subprocess.run([sys.executable, str(REPO / "tests" / script)],
                          cwd=str(REPO), env=env, capture_output=True, text=True)


@pytest.mark.parametrize("script", SELF_CONTAINED)
def test_standalone_regression(script):
    proc = _run(script)
    assert proc.returncode == 0, f"{script} failed:\n{proc.stdout}\n{proc.stderr}"


@pytest.mark.parametrize("script,required", sorted(NEEDS_KITEX_OUTPUT.items()))
def test_regression_against_example_output(script, required):
    if not (REPO / required).exists():
        pytest.skip(f"needs {required} -- run the generating script and KITEx first")
    proc = _run(script)
    assert proc.returncode == 0, f"{script} failed:\n{proc.stdout}\n{proc.stderr}"
