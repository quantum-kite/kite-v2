"""custom_two must run at the temperature the script asked for.

Python writes every `Temperature` dataset already divided by `energy_scale`,
and the readers' internal `temperature` is in those same rescaled units -- each
reader's own default (`0.001/scale`) and its own `-T` override (`.../scale`)
say so. `customtwo.cpp` used to divide the stored value a second time, so every
custom_two ran at `T/energy_scale` while `conductivity_dc` ran at `T` from the
identical number in the identical file.

The check is a comparison, not a stored constant: one file declares both
calculations at the same temperature, and KITE-tools has to report the same
value for both -- and that value has to be the one the script asked for. A
regression here is invisible in the output curve, which is why it survived; the
banner is the only place the two conventions ever meet.
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
FIXTURE = REPO / "tests" / "fixtures" / "graphene_custom_two_xx_check.py"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")

#: What the fixture asks for, and the scale it asks for it at.
TEMPERATURE = 0.01
ENERGY_SCALE = 4.0

_BANNER = re.compile(r"^\s*Temperature:\s*([\d.eE+-]+)(.*)$", re.MULTILINE)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("custom_two_temperature")
    subprocess.run([sys.executable, str(FIXTURE)], cwd=d, check=True,
                   capture_output=True, text=True)
    h5 = next(d.glob("*-output.h5"))
    subprocess.run([str(KITEX), h5.name], cwd=d, check=True, capture_output=True)
    return d, h5


def _banner(run, *extra):
    d, h5 = run
    out = subprocess.run([str(KITE_TOOLS), h5.name, "--CondDC", "-N", "cdc.dat", *extra],
                         cwd=d, check=True, capture_output=True, text=True).stdout
    return [(float(v), flag) for v, flag in _BANNER.findall(out)]


def test_the_fixture_still_states_the_values_this_test_assumes():
    """The comparison is only sharp while the two numbers differ by the scale."""
    source = FIXTURE.read_text()
    assert f"temperature={TEMPERATURE}" in source and f"temperature_={TEMPERATURE}" in source, (
        "the fixture no longer asks both calculations for the same temperature")
    assert f"spectrum_range=[-{ENERGY_SCALE:.0f}, {ENERGY_SCALE:.0f}]" in source, (
        "the fixture's energy scale changed, so a doubled division would no "
        "longer be distinguishable from the value asked for")


def test_both_calculations_store_the_same_rescaled_temperature(run):
    """The Python side writes one convention; the asymmetry was downstream."""
    _, h5 = run
    with h5py.File(h5, "r") as f:
        dc = float(np.array(f["Calculation/conductivity_dc/Temperature"]).ravel()[0])
        custom = float(np.array(f["Calculation/CustomTwo/Temperature"]).ravel()[0])
    assert dc == pytest.approx(custom)
    assert dc == pytest.approx(TEMPERATURE / ENERGY_SCALE)


def test_both_readers_report_the_temperature_that_was_asked_for(run):
    """The banner is where the two conventions meet, so it is where this is caught."""
    reported = _banner(run)
    assert len(reported) == 2, (
        f"expected one Temperature line per calculation, got {reported}")
    for value, flag in reported:
        assert "default" not in flag, (
            "a reader fell back to its default temperature instead of the "
            "stored one, so this run says nothing about the stored convention")
        assert value == pytest.approx(TEMPERATURE, rel=1e-6), (
            f"KITE-tools reports {value} where the script asked for "
            f"{TEMPERATURE}; a factor of {TEMPERATURE / value:.3g} means the "
            "stored value was rescaled one time too many or too few")


def test_the_command_line_override_agrees_with_the_stored_value(run):
    """`-T` already went through `/scale` in both readers; it must still match."""
    override = 0.02
    reported = _banner(run, "-T", str(override), "--CustomTwo", "-T", str(override))
    assert len(reported) == 2
    for value, _ in reported:
        assert value == pytest.approx(override, rel=1e-6)
