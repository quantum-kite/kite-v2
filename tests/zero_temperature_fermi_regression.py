"""A requested temperature of zero must give numbers, and must change nothing else.

`beta = 1/kT` is infinite when the user asks for T = 0, which is the Python
default for `conductivity_dc`, `conductivity_optical` and
`conductivity_optical_nonlinear`. Away from the chemical potential that is
already the right answer -- `exp` saturates and the occupation is 1 or 0 -- but
where an integration sample lands exactly on mu the product is `inf * 0`, which
is NaN, and one NaN poisons the whole Simpson sum rather than one point of it.

Two things are asserted here, and the second matters more than the first.

The repair is confined to inputs that previously produced NaN: `fermi_function`
returns 1/2 only when its exponent is NaN, and every other input takes the
expression it always took. So a run that produced numbers before must produce
the *same* numbers now. That was checked against the Haldane quantum-anomalous-
Hall and graphene quantum-Hall examples byte-for-byte before the change was
made, and the plateau test below is the part of that guarantee worth keeping in
CI: sigma_xy in the gap is quantised to the Chern number, and no change to the
occupation factor may disturb it.

KITE-tools parallelises its energy integral, and the reduction order is not
fixed, so repeated runs of the *same* binary differ in the last digits. Every
comparison here therefore pins `-t 1`; without it a byte comparison measures
thread scheduling, not the code under test.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")

#: Haldane with t2 = t/10 and the phase as the example codes it: C = +1.
CHERN = 1.0

_BUILDER = r'''
import sys, numpy as np, kite
sys.path.insert(0, sys.argv[2])
from dos_dccond_haldane import haldane   # the shipped model, not a copy of it

lat = haldane()
cfg = kite.Configuration(divisions=[2, 2], length=[128, 128],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=0, spectrum_range=[-10, 10],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
# num_points is the Fermi grid; odd on purpose, so it contains the band centre
# and a zero temperature puts a sample exactly on mu.
calc.conductivity_dc(num_points=201, num_moments=256, num_random=4,
                     num_disorder=1, direction='xy')   # temperature defaults to 0
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


#: The second-order reader is where the coincidence is structural rather than
#: incidental: its energy grid is LinSpaced(513, -0.995, 0.995) -- odd and
#: symmetric, so it contains exactly 0 -- and its default Fermi energy is
#: exactly 0 on an unshifted model. A CondDC grid only lands on mu for some
#: combinations of precision, range and point count, so it is not used here to
#: test the singular case.
_NONLINEAR = r'''
import sys, numpy as np, kite
from kite import lattice as latt

a, a_cc, t, gap = 0.246, 0.142, -2.8, 4.5
lat = latt.Lattice(a1=[a/2, a/2*np.sqrt(3)], a2=[a/2, -a/2*np.sqrt(3)])
lat.add_sublattices(('A', [0, 0], -gap/2.0), ('B', [0, a_cc], gap/2.0))
lat.add_hoppings(([0, 0], 'A', 'B', t), ([0, 1], 'A', 'B', t), ([-1, 0], 'A', 'B', t))

cfg = kite.Configuration(divisions=[2, 2], length=[64, 64],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-9.5, 9.5],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
# No temperature argument: the Python default is 0, which is what every shipped
# second-order script uses and what made beta infinite.
calc.conductivity_optical_nonlinear("xxy", num_points=64, num_moments=64,
                                    num_random=1, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def nonlinear(tmp_path_factory):
    d = tmp_path_factory.mktemp("zeroT_nlo")
    subprocess.run([sys.executable, "-c", _NONLINEAR, "n.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "n.h5"], cwd=d, check=True, capture_output=True)
    return d


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("zeroT")
    subprocess.run([sys.executable, "-c", _BUILDER, "h.h5", str(REPO / "examples")],
                   cwd=d, check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "h.h5"], cwd=d, check=True, capture_output=True)
    return d


def _conductivity(d, name, *extra):
    subprocess.run([str(KITE_TOOLS), "h.h5", "--CondDC", "-t", "1", "-N", name, *extra],
                   cwd=d, check=True, capture_output=True, text=True)
    return np.loadtxt(d / name)


def test_zero_temperature_produces_numbers(nonlinear):
    """The case where the sample lands on mu by construction, not by luck."""
    subprocess.run([str(KITE_TOOLS), "n.h5", "--CondOpt2", "-N", "nlo.dat"],
                   cwd=nonlinear, check=True, capture_output=True, text=True)
    data = np.loadtxt(nonlinear / "nlo.dat")
    bad = (~np.isfinite(data[:, 1:]).all(axis=1)).sum()
    assert bad == 0, (
        f"{bad} of {len(data)} rows are NaN at the default temperature of zero. "
        "The energy grid contains exactly 0 and so does the Fermi energy, so the "
        "occupation factor evaluates inf*0, and one NaN takes the whole Simpson "
        "sum with it rather than one point of it")


def test_the_quantum_anomalous_hall_plateau_is_quantised(run):
    """What the occupation factor must not be allowed to disturb."""
    data = _conductivity(run, "t0b.dat")
    energy, sigma_xy = data[:, 0], data[:, 1]
    gap = np.abs(energy) < 0.08
    assert gap.sum() > 10, "the window no longer selects points inside the gap"
    assert sigma_xy[gap].mean() == pytest.approx(CHERN, abs=0.05), (
        f"sigma_xy in the gap averages {sigma_xy[gap].mean():.4f} e^2/h, not the "
        f"Chern number {CHERN}")


def test_a_finite_temperature_result_is_untouched_by_the_zero_case(run):
    """The guard fires only on NaN, so finite-temperature runs are bit-identical.

    Asserted as reproducibility of the same binary at fixed thread count: this
    is the property that was compared across the change itself, and it is what
    would break if the guard were ever widened to catch finite exponents.
    """
    first = _conductivity(run, "warm1.dat", "-T", "0.05")
    second = _conductivity(run, "warm2.dat", "-T", "0.05")
    assert np.array_equal(first, second), (
        "a finite-temperature run is not reproducible at -t 1, so the "
        "byte-comparison this change was validated with no longer means anything")
    assert np.isfinite(first).all()
