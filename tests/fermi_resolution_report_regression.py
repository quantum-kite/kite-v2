"""A temperature the energy grid cannot resolve must be reported, not absorbed.

The Fermi factor is analytic with poles a distance pi*kT off the real axis, so
integrating it on a uniform grid converges geometrically in kT/dE, where dE is
the integration step. Measured on the shipped CondDC grid against a reference
integral: relative error 2.4e-7 at a ratio of 1, 1.8e-5 at 0.5, 1.6e-4 at 0.25,
and 1.3e-3 at zero temperature, where the integrand has a genuine step and
Simpson loses its order entirely.

Below roughly half a step, then, the temperature the user asked for is not the
one the calculation represents -- and the banner printed it back as if it were.
The readers now print the ratio, and say so when it is too small.

This is reported rather than enforced. A user may legitimately want a
temperature the grid cannot resolve; what they may not have is no way to tell.
The same argument is already made in conductivity_dc's own source comment for
the other length in the integrand, the broadening: "The step must stay well
under `scat`". This is that argument applied to the other factor.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")

_RATIO = re.compile(r"Fermi edge in integration steps:\s*([\d.eE+-]+)")
UNRESOLVED = "narrower than half the integration step"

_BUILDER = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

a, acc = 0.24595, 0.142
lat = latt.Lattice(a1=a*np.array([1.0, 0.0]), a2=a*np.array([0.5, 0.5*np.sqrt(3)]))
lat.add_sublattices(('A', [0, -acc/2], 0.0), ('B', [0, acc/2], 0.0))
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([1, -1], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0))

cfg = kite.Configuration(divisions=[2, 2], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-4, 4],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
calc.conductivity_dc(num_points=100, num_moments=64, num_random=1,
                     num_disorder=1, direction='xx', temperature=0.01)
A = custom.Vertex(64, [[1.0, "vx"]])
B = custom.Vertex(64, [[1.0, "vx"]])
calc.custom_two(stream_=[A, B], num_random_=1, num_disorder_=1,
                num_points_=100, temperature_=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("fermi_res")
    subprocess.run([sys.executable, "-c", _BUILDER, "r.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "r.h5"], cwd=d, check=True, capture_output=True)

    def tools(*extra):
        return subprocess.run([str(KITE_TOOLS), "r.h5", "--CondDC", "-t", "1",
                               "-N", "o.dat", *extra],
                              cwd=d, check=True, capture_output=True, text=True).stdout
    return d, tools


def test_the_ratio_is_reported(run):
    """Every reader with a Fermi factor and an energy grid must print it."""
    _, tools = run
    found = _RATIO.findall(tools())
    assert len(found) == 2, (
        f"expected the ratio from both conductivity_dc and custom_two, got {found}")


def test_a_resolvable_temperature_is_not_flagged(run):
    """The marker must stay quiet when the grid is adequate, or it is noise."""
    _, tools = run
    out = tools("-T", "0.5")
    ratios = [float(x) for x in _RATIO.findall(out)]
    assert all(r > 0.5 for r in ratios), ratios
    assert UNRESOLVED not in out


def test_zero_temperature_is_flagged(run):
    """The case the readers used to accept silently."""
    _, tools = run
    out = tools("-T", "0")
    assert [float(x) for x in _RATIO.findall(out)][0] == pytest.approx(0.0)
    assert UNRESOLVED in out, (
        "a temperature of zero no longer reports that the Fermi edge is "
        "unresolved, so the slow convergence is invisible again")


def test_a_small_but_nonzero_temperature_is_flagged(run):
    """Not only the singular case: anything under half a step is unresolved."""
    _, tools = run
    # dE = 1.98/1024 in reduced units, so this is well under half a step.
    out = tools("-T", "0.0005")
    assert [float(x) for x in _RATIO.findall(out)][0] < 0.5
    assert UNRESOLVED in out


def test_the_report_does_not_change_the_numbers(run):
    """It is a diagnostic. It must not touch the result it describes."""
    d, tools = run
    tools()
    first = (d / "o.dat").read_bytes()
    tools()
    assert (d / "o.dat").read_bytes() == first
