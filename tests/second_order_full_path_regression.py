"""The general second-order path must run, and the example must ask for it.

`conductivity_optical_nonlinear` has two code paths. The general one computes
four Chebyshev objects, Gamma0 through Gamma3; an internal shortcut computes
only Gamma1 and Gamma2 and is a materials-specific simplification, not a
numerical one. Until this cycle the only shipped script exercising the
second-order code requested the shortcut, so the general path had no example
and no test, and the one result file in the tree was NaN throughout.

The NaN was not the path: it was a requested temperature of zero making
`beta` infinite (see tests/zero_temperature_fermi_regression.py). With that
fixed the general path runs, and this pins it.

What is asserted is reproducibility and structure, not a stored curve. A
second-order conductivity at a size CI can afford is nowhere near converged --
the peak moves by a factor of two between `num_random` 4 and 16 -- so a golden
file would pin noise. Two independent exports with the same seeds must agree
bit-for-bit, all four Gammas must be present, and the response must be finite.
"""

import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"
FIXTURE = REPO / "tests" / "fixtures" / "hbn_optcond2_full_path_check.py"
EXAMPLE = REPO / "examples" / "hbn_optcond2_vacancies.py"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")

GROUP = "Calculation/conductivity_optical_nonlinear"


def _build(d, name):
    subprocess.run([sys.executable, str(FIXTURE), name], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), name], cwd=d, check=True, capture_output=True)
    return d / name


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("nlo_full")
    h5 = _build(d, "a.h5")
    subprocess.run([str(KITE_TOOLS), h5.name, "--CondOpt2", "-N", "a.dat"],
                   cwd=d, check=True, capture_output=True, text=True)
    return d, h5


def test_all_four_chebyshev_objects_are_computed(run):
    """The shortcut writes only Gamma1 and Gamma2; the general path writes four."""
    _, h5 = run
    with h5py.File(h5, "r") as f:
        present = sorted(k for k in f[GROUP] if k.startswith("Gamma"))
    assert len(present) == 4, (
        f"expected Gamma0..Gamma3, found {present}. Fewer than four means the "
        "run took the two-object shortcut, so the general path is untested")


def test_the_response_is_finite(run):
    """The acceptance criterion: the full path runs end to end."""
    d, _ = run
    data = np.loadtxt(d / "a.dat")
    bad = (~np.isfinite(data[:, 1:]).all(axis=1)).sum()
    assert bad == 0, f"{bad} of {len(data)} rows are not finite"
    assert np.abs(data[:, 1:]).max() > 0, (
        "the whole response is identically zero, which for a gapped honeycomb "
        "means the inversion-breaking onsite term is not reaching the vertex")


def test_two_exports_with_the_same_seeds_agree_exactly(run):
    """A fixture is only worth storing if the thing it stores is deterministic."""
    d, _ = run
    second = _build(d, "b.h5")
    subprocess.run([str(KITE_TOOLS), second.name, "--CondOpt2", "-N", "b.dat"],
                   cwd=d, check=True, capture_output=True, text=True)
    assert (d / "a.dat").read_bytes() == (d / "b.dat").read_bytes(), (
        "two exports with identical seed_h/seed_v produced different results, "
        "so the second-order path is not reproducible and cannot be pinned")


def test_the_shipped_example_requests_the_general_path():
    """AC-NLO-01's precondition: the example must exercise what is claimed."""
    source = EXAMPLE.read_text()
    call = source[source.index("conductivity_optical_nonlinear("):]
    call = call[:call.index(")")]
    assert "special" not in call, (
        "hbn_optcond2_vacancies.py asks for the two-object shortcut again, so "
        "the general second-order path has no shipped example")
    assert re.search(r'"xxy"', call), (
        "the example no longer requests the xxy component its docstring explains")
