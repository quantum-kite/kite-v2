"""The default integration grid must resolve the default broadening.

Both Kubo-Bastin reconstructions integrate over energy on a fixed grid, and
both default their broadening to 0.01 eV. The grid step has to sit well under
that broadening or the Lorentzian is under-sampled and the answer is wrong by a
percent or so -- quietly, since nothing about the output looks unusual.

At the former default of 513 points over [-0.99, 0.99] the step was 3.87e-3
against a broadening of 0.01/4.0 = 2.5e-3: coarser than the feature being
integrated. Measured against a converged grid:

    CustomTwo, Kane-Mele plateau   513 -> -2.0656   2049 -> -2.0874   (1.05% off)
    CondDC, graphene sigma_xx      513 ->  6.9859   2049 ->  6.8521   (1.96% off)

so the default is now 2049, which agrees with a 32769-point grid to 2e-4%.

The arithmetic check below is the cheap half and runs anywhere. The end-to-end
half needs a built KITE-tools and stored moments, so it skips without them.
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITE_TOOLS = REPO / "build" / "KITE-tools"
SWEEP = REPO / "sweep-data"

#: Each reconstruction's default grid, read from its own source.
SOURCES = {
    "custom_two": REPO / "tools" / "Src" / "Custom" / "customtwo.cpp",
    "conductivity_dc": REPO / "tools" / "Src" / "CondDC" / "conductivity_dc.cpp",
}

#: The integration window both reconstructions use, and the default broadening
#: in KPM reduced units for the examples here (0.01 eV over an energy_scale of 4).
WINDOW = 1.98
BROADENING = 0.01 / 4.0


def _default_grid(path):
    match = re.search(r"NEnergies\s*=\s*(\d+)", path.read_text())
    assert match, f"{path.name} no longer sets a default NEnergies"
    return int(match.group(1))


@pytest.mark.parametrize("name", sorted(SOURCES))
def test_the_default_step_is_finer_than_the_default_broadening(name):
    """dE << eta, or the integrator samples a Lorentzian it cannot see."""
    points = _default_grid(SOURCES[name])
    step = WINDOW / (points - 1)
    assert step < BROADENING / 2, (
        f"{name}: the default grid of {points} points gives dE = {step:.2e}, "
        f"against a default broadening of {BROADENING:.2e}. The step must stay "
        f"well under the broadening or the energy integral is under-resolved.")


@pytest.mark.skipif(not KITE_TOOLS.exists(), reason="needs a built KITE-tools")
@pytest.mark.parametrize("example,flag,output,expected", [
    ("kane_mele_spin_hall", "--CustomTwo", "custom_two.dat", -2.087378),
    ("graphene_custom_two_xx_check", "--CondDC", "condDC.dat", 6.851867),
])
def test_the_default_reproduces_the_converged_result(example, flag, output, expected):
    """End to end: running with no -E must match a 32769-point grid.

    This is what the arithmetic above is a proxy for, and it is the assertion
    that would have failed before the default was raised.
    """
    directory = SWEEP / example
    h5 = next(directory.glob("*.h5"), None)
    if h5 is None:
        pytest.skip(f"needs {directory.name}/; run tests/run_all_examples.py first")

    subprocess.run([str(KITE_TOOLS), h5.name, flag], cwd=str(directory),
                   capture_output=True, text=True, check=True)
    data = np.loadtxt(directory / output)
    plateau = data[np.abs(data[:, 0]) < 0.6, 1].mean()
    assert plateau == pytest.approx(expected, rel=1e-3), (
        f"the default grid gives {plateau:.6f} where a converged grid gives "
        f"{expected:.6f}; the default may have been lowered again")
