"""The analytic moment ladder: a correctness check that needs no physics input.

The first few Chebyshev moments of the rescaled Hamiltonian are computable from
the lattice specification alone, with no reference to what KITE produced:

    mu_0 = 1                                  (the DOS sum rule)
    mu_1 = (mean onsite - shift) / scale
    mu_2 = 2 * Tr(H~^2)/N - 1
    mu_odd = 0   on any bipartite lattice with zero onsite energy (chiral symmetry)

Comparing KITE's `/Calculation/dos/MU` against those values tests, in one pass
and for any model: the sum rule, the hopping magnitudes, the energy scale and
shift, the disorder statistics, and chiral symmetry. None of it requires
knowing what the calculation is *for*.

It does **not** test the Peierls phases, and an earlier version of this
docstring claimed it did. `mu_2` depends on `Tr H^2 = sum |H_ij|^2`, which is
blind to the phase of every hopping, and the odd moments vanish on a bipartite
lattice by chiral symmetry whatever the flux. The first moment that sees a flux
is `mu_3`, which closes a triangle, and the ladder stops at `mu_2`. The phases
are covered by `magnetic_field_unit_cell_area_regression.py`.

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

import os
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


#: The lattice the runs below use, in unit cells. The site count follows from
#: this and the sublattice count, which is why nothing here hardcodes it.
LENGTH = [64, 64]


def shifted_square(t=1.0, onsite=1.5):
    """A square lattice with a non-zero onsite term, so `energy_shift != 0`.

    Both lattices above have zero onsite energy and are run on ranges symmetric
    about zero, which makes `energy_shift` exactly 0 in every case: `mu_1` is
    then identically zero, the `-2*shift*mean_onsite + shift**2` terms in
    `analytic_moments` are never evaluated, and a sign error on the shift or a
    shift/scale mix-up is invisible. This lattice exists to make them visible.
    """
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(("A", [0, 0], onsite))
    lat.add_hoppings(([1, 0], "A", "A", -t), ([0, 1], "A", "A", -t))
    return lat


def run_moments(lattice, tmp_path, spectrum_range, num_random=8, num_moments=32,
                disorder=None):
    """Export, run KITEx, and return the moments it wrote."""
    configuration = kite.Configuration(
        divisions=[1, 1], length=LENGTH, boundaries=["periodic", "periodic"],
        is_complex=False, precision=1, spectrum_range=spectrum_range,
        seed_h=20261101, seed_v=20261101)
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=200, num_moments=num_moments,
                    num_random=num_random, num_disorder=1)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "ladder.h5"
    kwargs = {"disorder": disorder} if disorder is not None else {}
    kite.config_system(lattice, configuration, calculation, filename=str(path), **kwargs)
    # Not SEED=ones: a vector of ones is a particular vector, not a random one,
    # and the stochastic trace stops estimating Tr/N -- mu_1 comes back at -0.91
    # instead of 0. The probe vectors are pinned through `seed_v` instead, which
    # keeps them random and keeps two runs drawing the same ones.
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    subprocess.run([str(KITEX), path.name], cwd=tmp_path, check=True,
                   capture_output=True, env=environment)
    with h5py.File(path, "r") as f:
        return (np.asarray(f["Calculation/dos/MU"][()]).ravel(),
                float(f["EnergyScale"][()]), float(f["EnergyShift"][()]))


def sites(lattice):
    """Sites in the run, counting the orbital index -- cells times sublattices."""
    return int(np.prod(LENGTH)) * lattice.nsub


def tolerance(lattice, num_random, k=5.0):
    """5 / sqrt(R*N): the stochastic-trace floor, not a fudge factor.

    The floor is only a floor if N is the site count of the run being judged.
    Taking it from the lattice rather than writing a number keeps that true:
    graphene has two sublattices and the square lattice one, so a fixed N would
    hold one of them to a tolerance sqrt(2) away from its own floor and the
    docstring above would stop being accurate.
    """
    return k / np.sqrt(num_random * sites(lattice))


# --------------------------------------------------------------------------
# the tests
# --------------------------------------------------------------------------

@pytest.mark.parametrize("build,half_width", [(graphene, 3.3), (square, 4.4)])
def test_mu0_is_one(build, half_width, tmp_path):
    """The DOS sum rule, per site including the orbital index."""
    lattice = build()
    mu, _, _ = run_moments(lattice, tmp_path, [-half_width, half_width])
    assert mu[0] == pytest.approx(1.0, abs=tolerance(lattice, 8))


@pytest.mark.parametrize("build,half_width", [(graphene, 3.3), (square, 4.4)])
def test_mu1_and_mu2_match_the_lattice(build, half_width, tmp_path):
    lattice = build()
    mu, scale, shift = run_moments(lattice, tmp_path, [-half_width, half_width])
    _, mu1, mu2 = analytic_moments(lattice, scale, shift)
    assert mu[1] == pytest.approx(mu1, abs=tolerance(lattice, 8))
    assert mu[2] == pytest.approx(mu2, abs=tolerance(lattice, 8))


def test_odd_moments_vanish_on_a_bipartite_lattice(tmp_path):
    """Exact by chiral symmetry, so any violation is unambiguous."""
    lattice = graphene()
    mu, _, _ = run_moments(lattice, tmp_path, [-3.3, 3.3])
    odd = mu[1::2]
    assert np.max(np.abs(odd)) < tolerance(lattice, 8), (
        f"largest odd moment {np.max(np.abs(odd)):.2e} exceeds the stochastic floor")


def test_disorder_enters_mu2_as_its_variance(tmp_path):
    """Tr(H^2)/N gains exactly sigma^2 per site.

    Asserted as a *difference* between a clean run and a disordered one on the
    same lattice, at the same probe vector, rather than as an absolute value.
    The whole effect is `2*sigma^2/scale^2 = 0.0113`, which is a fifth of the
    stochastic floor this file tolerates elsewhere: an absolute comparison here
    passes just as happily if KITE drops the disorder term altogether, which is
    the one thing the test exists to notice.

    What the difference distinguishes is the three readings of `sigma` that
    matter -- the standard deviation, a box half-width (which would give
    sigma^2/3), and no disorder at all. Rather than invent a tolerance, the
    assertion is that the measured difference is closest to the first.
    """
    lattice = graphene()
    sigma = 0.3
    disorder = kite.Disorder(lattice)
    disorder.add_disorder("A", "Uniform", 0.0, sigma)
    disorder.add_disorder("B", "Uniform", 0.0, sigma)

    clean, scale, _ = run_moments(lattice, tmp_path / "clean", [-4.0, 4.0],
                                  num_random=32)
    dirty, _, _ = run_moments(lattice, tmp_path / "dirty", [-4.0, 4.0],
                              num_random=32, disorder=disorder)
    measured = dirty[2] - clean[2]

    candidates = {
        "sigma is the standard deviation": 2 * sigma ** 2 / scale ** 2,
        "sigma is a box half-width": 2 * (sigma ** 2 / 3) / scale ** 2,
        "the disorder never reached the Hamiltonian": 0.0,
    }
    best = min(candidates, key=lambda k: abs(measured - candidates[k]))
    assert best == "sigma is the standard deviation", (
        f"mu_2 moved by {measured:.5f} between the clean and disordered runs; "
        f"that is closest to '{best}' ({candidates[best]:.5f}), not to "
        f"{candidates['sigma is the standard deviation']:.5f}")
    assert measured == pytest.approx(
        candidates["sigma is the standard deviation"], rel=0.25), (
        f"the shift in mu_2 is {measured:.5f} against an expected "
        f"{candidates['sigma is the standard deviation']:.5f}")


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


def test_the_shift_terms_are_exercised_at_all(tmp_path):
    """An asymmetric range on a lattice with an onsite term, so shift != 0."""
    lattice = shifted_square(onsite=1.5)
    # The spectrum is [-2.5, 5.5]; the range is wider and deliberately NOT
    # centred on it, and its centre 2.5 is deliberately NOT the onsite term --
    # picking the two equal would put mu_1 back at zero and lose the point.
    mu, scale, shift = run_moments(lattice, tmp_path, [-3.0, 8.0])
    assert abs(shift) > 0.1, (
        f"energy_shift came back as {shift}; this case exists to make the "
        "shift terms non-trivial and it is not doing so")
    assert abs(mu[1]) > 0.1, (
        f"mu_1 is {mu[1]}, still effectively zero; the first moment is supposed "
        "to carry (mean onsite - shift)/scale here")
    _, mu1, mu2 = analytic_moments(lattice, scale, shift)
    assert mu[1] == pytest.approx(mu1, abs=tolerance(lattice, 8))
    assert mu[2] == pytest.approx(mu2, abs=tolerance(lattice, 8))
