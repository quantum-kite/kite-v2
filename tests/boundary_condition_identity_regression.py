"""Two boundary conditions KITE can check exactly, and never did.

A twist enters as a phase on the bond that wraps the sample
(`Src/Vector/KPM_Vector2D.cpp:195-205`: `exp(i*phi)` guarded by an integer
division that is 1 only for the wrapping hop). Two consequences are exact, not
approximate, and both were available and unused:

* a twist of zero must reproduce periodic boundaries **bit for bit**, since the
  phase is then 1 and the Hamiltonian is the same matrix;
* a twist of pi multiplies the wrapping bond by -1, which *is* antiperiodic.

The second needs a reference built without the twist machinery, or it only
compares the code to itself. Antiperiodic boundaries on a ring of `L` cells are
related to periodic ones by the gauge `psi_n -> exp(i*pi*n_x/L) psi_n`, which
moves the sign off the wrapping bond and spreads it as `exp(i*pi/L)` over every
bond along that axis. So the reference here is an ordinary **periodic** lattice
with a complex hopping, and no twist at all.

The observable is the local density of states rather than the total, because the
comparison has to survive that gauge. A diagonal gauge leaves `<j|T_n(H)|j>`
alone -- the two phases cancel -- while a stochastic trace over random vectors
does not, since the probe vector would have to be gauge-transformed too. The LDoS
is also seeded from a basis ket, so with no disorder there is nothing random left
and "to machine precision" means what it says.

`open` and `random` are exercised too, but only the two identities above are
exact. A random twist is drawn per vector and is not reproducible from the inputs
alone, which is its own subject.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

L = 16
PROBE = [[8, 8]]

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

L = 16
name, boundaries, angles, tx = (
    sys.argv[1], sys.argv[2].split(','), [float(a) for a in sys.argv[3].split(',')],
    complex(sys.argv[4]))

lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', tx), ([0, 1], 'A', 'A', -1.0))

cfg = kite.Configuration(divisions=[1, 1], length=[L, L], boundaries=boundaries,
                         angles=angles, is_complex=True, precision=1,
                         spectrum_range=[-5, 5])
calc = kite.Calculation(cfg)
calc.ldos(energy=[0.0], num_moments=64, num_disorder=1,
          position=[[8, 8]], sublattice=['A'])
kite.config_system(lat, cfg, calc, filename=name)
'''

#: name -> (boundaries, angles, hopping along a1)
CASES = {
    "periodic":     ("periodic,periodic", "0,0", "-1"),
    "twist_zero":   ("twisted,periodic",  "0,0", "-1"),
    "twist_pi":     ("twisted,periodic",  f"{np.pi},0", "-1"),
    # antiperiodic, written in the gauge that needs no twist at all
    "antiperiodic": ("periodic,periodic", "0,0", str(-np.exp(1j * np.pi / L))),
    "open":         ("open,periodic",     "0,0", "-1"),
    "random":       ("random,periodic",   "0,0", "-1"),
}


@pytest.fixture(scope="module")
def moments(tmp_path_factory):
    import h5py
    d = tmp_path_factory.mktemp("boundaries")
    out = {}
    for name, (boundaries, angles, tx) in CASES.items():
        h5 = d / f"{name}.h5"
        subprocess.run([sys.executable, "-c", _BUILDER, str(h5), boundaries,
                        angles, tx], cwd=d, check=True, capture_output=True,
                       text=True)
        subprocess.run([str(KITEX), str(h5)], cwd=d, check=True,
                       capture_output=True)
        with h5py.File(h5, "r") as f:
            out[name] = np.asarray(f["/Calculation/ldos/lMU"][()]).ravel()
    return out


def test_a_zero_twist_is_periodic_bit_for_bit(moments):
    """The phase is 1, so this is the same matrix, not merely a close one."""
    diff = np.abs(moments["twist_zero"] - moments["periodic"]).max()
    assert diff == 0.0, (
        f"a twist of zero differs from periodic by {diff:.3e}. The twist enters "
        "as exp(i*phi) on the wrapping bond, so at phi = 0 the Hamiltonian is "
        "identical and the moments must be too")


def test_a_twist_of_pi_is_antiperiodic(moments):
    """Against a reference built with a complex hopping and no twist."""
    diff = np.abs(moments["twist_pi"] - moments["antiperiodic"]).max()
    assert diff < 1e-13, (
        f"a twist of pi differs from antiperiodic by {diff:.3e}. The twist puts "
        "-1 on the wrapping bond; the reference spreads the same sign over every "
        "bond along that axis as exp(i*pi/L), which is the same Hamiltonian in a "
        "diagonal gauge the LDoS cannot see")


def test_the_identities_are_not_vacuous(moments):
    """A twist of pi must actually change the answer, or the above proves nothing."""
    spread = np.abs(moments["twist_pi"] - moments["periodic"]).max()
    assert spread > 1e-3, (
        f"periodic and twist-pi agree to {spread:.3e}, so this system cannot "
        "tell the two boundary conditions apart and neither identity above "
        "carries any weight")


def test_open_boundaries_are_distinct_and_finite(moments):
    """The fourth boundary mode; cutting the wrapping bond is not a phase."""
    assert np.isfinite(moments["open"]).all()
    assert np.abs(moments["open"] - moments["periodic"]).max() > 1e-3


def test_random_twists_run_and_stay_finite(moments):
    """Smoke only: the draw is not reproducible from the inputs alone."""
    assert np.isfinite(moments["random"]).all()
    assert moments["random"][0] == pytest.approx(1.0), (
        "the zeroth moment is 1 by normalisation whatever the boundary is")
