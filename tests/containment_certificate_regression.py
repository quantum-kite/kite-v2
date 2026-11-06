"""KITEx must stop when the rescaled spectrum is not inside [-1, 1].

Chebyshev polynomials are bounded by 1 on [-1, 1], so for a spectrum that
really is contained, `||phi_n|| <= ||phi_0||` at every step. If it escapes, the
recursion grows exponentially -- and KITEx used to finish anyway, exit 0, and
write moments of 1e85 against a physical bound of 1. The number is meaningless
long before the run ends, and nothing said so.

`spectrum_range` cannot be made rigorous at export: with Gaussian
disorder the support is unbounded, so there is no range that strictly contains
the spectrum. The certificate is the other half -- it costs one norm every 16
steps and it is exact, because the bound holds if and only if the spectrum is
contained.

The threshold is 1.01 and not 1: an adequate range sits at exactly 1.000, so
there is room to spare and no false positive. The whole test suite running
clean is the real evidence for that, since every other example here has a range
the author chose.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

# Clean square lattice: the band is exactly [-4, 4], no disorder, so whether
# the range contains the spectrum is a fact and not an estimate.
lo, hi = float(sys.argv[2]), float(sys.argv[3])
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', -1.0), ([0, 1], 'A', 'A', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[lo, hi])
calc = kite.Calculation(cfg)
calc.dos(num_points=512, num_moments=256, num_random=1, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _run(tmp, lo, hi):
    path = tmp / f"r{lo}_{hi}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path), str(lo), str(hi)],
                   cwd=tmp, check=True, capture_output=True, text=True)
    return path, subprocess.run([str(KITEX), str(path)], capture_output=True,
                                text=True)


@pytest.fixture(scope="module")
def too_small(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("uncontained"), -2, 2)


@pytest.fixture(scope="module")
def adequate(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("contained"), -4.5, 4.5)


def test_an_uncontained_spectrum_aborts(too_small):
    _, proc = too_small
    assert proc.returncode == 1, (
        "a range of [-2, 2] cannot hold a band of [-4, 4]; the run must stop "
        f"rather than produce numbers (exit {proc.returncode})")


def test_the_abort_says_what_is_wrong_and_what_to_change(too_small):
    _, proc = too_small
    out = proc.stdout + proc.stderr
    assert "spectrum_range" in out, (
        f"the message does not name the parameter to change: {out.strip()!r}")
    assert "contained" in out.lower()


def test_it_stops_before_writing_the_moments(too_small):
    """The point is to fail before producing numbers, not after."""
    h5py = pytest.importorskip("h5py")
    path, _ = too_small
    with h5py.File(path, "r") as f:
        assert "/Calculation/dos/MU" not in f, (
            "the run wrote its moments anyway; aborting after the fact is no "
            "better than the exit-0 it replaced")


def test_an_adequate_range_is_untouched(adequate):
    """1.01 leaves room: a contained spectrum sits at exactly 1.000."""
    path, proc = adequate
    assert proc.returncode == 0, (
        f"a range of [-4.5, 4.5] contains the band and must run: {proc.stdout}")
    h5py = pytest.importorskip("h5py")
    import numpy as np
    with h5py.File(path, "r") as f:
        mu = np.asarray(f["/Calculation/dos/MU"][()]).ravel()
    # Not a precision check -- one random vector leaves several percent of
    # stochastic noise on mu_0. What matters is that the moments are O(1)
    # rather than the 1e85 the uncontained case used to write out.
    assert np.isfinite(mu).all(), "the moments are not finite"
    assert np.abs(mu).max() < 10.0, (
        f"the moments reach {np.abs(mu).max():.3g}; Chebyshev moments of a "
        "contained spectrum are bounded by 1")
