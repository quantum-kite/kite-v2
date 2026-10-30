"""The analytic moment ladder: a correctness check that needs no physics input.

The first few Chebyshev moments of the rescaled Hamiltonian are computable from
the lattice specification alone, with no reference to what KITE produced:

    mu_0 = 1                                  (the DOS sum rule)
    mu_1 = (mean onsite - shift) / scale
    mu_2 = 2 * Tr(H~^2)/N - 1
    mu_odd = 0   on any bipartite lattice with zero onsite energy (chiral symmetry)

Comparing KITE's `/Calculation/dos/MU` against those values tests, in one pass
and for any model: the sum rule, that the hoppings and Peierls phases were
written correctly, the energy scale and shift, the disorder statistics, and
chiral symmetry. None of it requires knowing what the calculation is *for*.

This is the check that would have caught F-63 on its own. When
`ldos_graphene.py` declared a spectral range a third of the true bandwidth, the
rescaled spectrum reached |x| = 3 and the odd moments came out around 1e48
instead of zero -- no knowledge of vacancies or LDOS needed.

Where it does NOT hold, and the tests say so rather than leaving it implicit:

  * the `green` (Lorentzian) kernel has algebraic tails and no sum rule;
  * with `is_complex=False` the sum rule holds only in expectation, with
    sigma(mu_0) = sqrt(0.8 / (N_eff * R));
  * peak height and width are set by the kernel, not the physics, and are never
    acceptance criteria.
"""

import subprocess

import h5py
import numpy as np
import pytest

import kite
from kite import lattice as latt

REPO = __import__("pathlib").Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
A_CC = 1 / np.sqrt(3)

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="KITEx not built")


# --------------------------------------------------------------------------
# the analytic side -- computed from the Lattice, never from KITE's output
# --------------------------------------------------------------------------

def analytic_moments(lattice, energy_scale, energy_shift=0.0, disorder_variance=0.0):
    """Return (mu_0, mu_1, mu_2) for the rescaled Hamiltonian.

    `disorder_variance` is the per-site variance added by on-site disorder --
    sigma^2 for a distribution of standard deviation sigma, since Tr(H^2)/N
    gains exactly the variance.
    """
    nsub = lattice.nsub
    onsite = np.zeros(nsub)
    for sub in lattice.sublattices.values():
        onsite[sub.alias_id] = float(np.real(np.asarray(sub.energy).ravel()[0]))

    # sum of |t|^2 over the bonds touching each sublattice; each stored term
    # also implies its conjugate, so it contributes to both endpoints
    row_sq = np.zeros(nsub)
    for family in lattice.hoppings.values():
        magnitude = float(np.abs(np.asarray(family.energy).ravel()[0])) ** 2
        for term in family.terms:
            row_sq[term.from_id] += magnitude
            row_sq[term.to_id] += magnitude

    mean_onsite = onsite.mean()
    trace_h2 = np.mean(onsite ** 2 + row_sq) + disorder_variance

    mu1 = (mean_onsite - energy_shift) / energy_scale
    h2_scaled = (trace_h2 - 2 * energy_shift * mean_onsite + energy_shift ** 2) / energy_scale ** 2
    return 1.0, mu1, 2 * h2_scaled - 1


def graphene(t=1.0):
    lat = latt.Lattice(a1=[1, 0], a2=[0.5, np.sqrt(3) / 2])
    lat.add_sublattices(("A", [0, -A_CC / 2], 0.0), ("B", [0, A_CC / 2], 0.0))
    lat.add_hoppings(([0, 0], "A", "B", -t), ([-1, 0], "A", "B", -t),
                     ([-1, 1], "A", "B", -t))
    return lat


def square(t=1.0):
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(("A", [0, 0], 0.0))
    lat.add_hoppings(([1, 0], "A", "A", -t), ([0, 1], "A", "A", -t))
    return lat


def run_moments(lattice, tmp_path, spectrum_range, num_random=8, num_moments=32,
                disorder=None):
    """Export, run KITEx, and return the moments it wrote."""
    configuration = kite.Configuration(
        divisions=[1, 1], length=[64, 64], boundaries=["periodic", "periodic"],
        is_complex=False, precision=1, spectrum_range=spectrum_range)
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=200, num_moments=num_moments,
                    num_random=num_random, num_disorder=1)
    path = tmp_path / "ladder.h5"
    kwargs = {"disorder": disorder} if disorder is not None else {}
    kite.config_system(lattice, configuration, calculation, filename=str(path), **kwargs)
    subprocess.run([str(KITEX), path.name], cwd=tmp_path, check=True, capture_output=True)
    with h5py.File(path, "r") as f:
        return (np.asarray(f["Calculation/dos/MU"][()]).ravel(),
                float(f["EnergyScale"][()]), float(f["EnergyShift"][()]))


def tolerance(num_random, num_sites=2 * 64 * 64, k=5.0):
    """5 / sqrt(R*N): the stochastic-trace floor, not a fudge factor."""
    return k / np.sqrt(num_random * num_sites)


# --------------------------------------------------------------------------
# the tests
# --------------------------------------------------------------------------

@pytest.mark.parametrize("build,half_width", [(graphene, 3.3), (square, 4.4)])
def test_mu0_is_one(build, half_width, tmp_path):
    """The DOS sum rule, per site including the orbital index."""
    mu, _, _ = run_moments(build(), tmp_path, [-half_width, half_width])
    assert mu[0] == pytest.approx(1.0, abs=tolerance(8))


@pytest.mark.parametrize("build,half_width", [(graphene, 3.3), (square, 4.4)])
def test_mu1_and_mu2_match_the_lattice(build, half_width, tmp_path):
    lattice = build()
    mu, scale, shift = run_moments(lattice, tmp_path, [-half_width, half_width])
    _, mu1, mu2 = analytic_moments(lattice, scale, shift)
    assert mu[1] == pytest.approx(mu1, abs=tolerance(8))
    assert mu[2] == pytest.approx(mu2, abs=tolerance(8))


def test_odd_moments_vanish_on_a_bipartite_lattice(tmp_path):
    """Exact by chiral symmetry, so any violation is unambiguous."""
    mu, _, _ = run_moments(graphene(), tmp_path, [-3.3, 3.3])
    odd = mu[1::2]
    assert np.max(np.abs(odd)) < tolerance(8), (
        f"largest odd moment {np.max(np.abs(odd)):.2e} exceeds the stochastic floor")


def test_disorder_enters_mu2_as_its_variance(tmp_path):
    """Tr(H^2)/N gains exactly sigma^2 per site."""
    lattice = graphene()
    disorder = kite.Disorder(lattice)
    sigma = 0.3
    disorder.add_disorder("A", "Uniform", 0.0, sigma)
    disorder.add_disorder("B", "Uniform", 0.0, sigma)
    mu, scale, shift = run_moments(lattice, tmp_path, [-4.0, 4.0], disorder=disorder)
    _, _, mu2 = analytic_moments(lattice, scale, shift, disorder_variance=sigma ** 2)
    assert mu[2] == pytest.approx(mu2, abs=3 * tolerance(8))


def test_the_ladder_catches_a_spectrum_range_that_is_too_narrow(tmp_path):
    """The F-63 failure -- now caught by KITEx before any moment is written.

    A range covering a third of the true bandwidth puts the rescaled spectrum
    at |x| = 3, where T_n grows exponentially. This test used to read the
    moments back and show they had exploded and that the odd ones no longer
    vanished. Since T-024 the recursion carries a containment certificate and
    the run aborts instead, so there are no moments to read: the same defect,
    detected earlier. See containment_certificate_regression.py.
    """
    with pytest.raises(subprocess.CalledProcessError):
        run_moments(graphene(t=1.0), tmp_path, [-1.0, 1.0], num_moments=24)
