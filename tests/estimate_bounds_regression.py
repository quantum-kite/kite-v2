"""`kite.estimate_bounds` must exist, work without pybinding, and be a bound.

It was documented as public API in `docs/api/kite.md` while never being imported
into the package, so calling it raised `AttributeError`. The old implementation
also built a pybinding model and read that library's KPM scaling factors, so it
could not work at all once pybinding became optional.

The replacement uses the Bloch Hamiltonian the lattice already defines, so the
clean-lattice term is *exact* rather than estimated, and adds the disorder
contributions additively -- legitimate by Weyl's inequality, and deliberately
over-estimating, which is the safe direction.

What matters is the direction of the error: the returned interval must CONTAIN
the spectrum. A bound that is slightly too wide costs resolution; one that is
too narrow puts eigenvalues outside [-1, 1], where the Chebyshev recursion
diverges and the output is meaningless with no warning.
"""

import numpy as np
import pytest

import kite
from kite import lattice as latt

A_CC = 1 / np.sqrt(3)


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


def cubic(t=1.0):
    lat = latt.Lattice(a1=[1, 0, 0], a2=[0, 1, 0], a3=[0, 0, 1])
    lat.add_sublattices(("A", [0, 0, 0], 0.0))
    lat.add_hoppings(([1, 0, 0], "A", "A", -t), ([0, 1, 0], "A", "A", -t),
                     ([0, 0, 1], "A", "A", -t))
    return lat


def test_it_is_importable_from_the_package():
    """The defect itself: documented as kite.estimate_bounds, never exported."""
    assert callable(kite.estimate_bounds)


def test_it_does_not_need_pybinding():
    import sys
    kite.estimate_bounds(graphene())
    assert "pybinding" not in sys.modules


@pytest.mark.parametrize("build,half_width", [
    (graphene, 3.0),        # +-3t at Gamma
    (square, 4.0),          # +-4t
    (cubic, 6.0),           # +-6t
])
def test_clean_band_extrema_are_exact(build, half_width):
    lo, hi = kite.clean_band_extrema(build())
    assert lo == pytest.approx(-half_width, abs=1e-6)
    assert hi == pytest.approx(+half_width, abs=1e-6)


@pytest.mark.parametrize("build,half_width", [(graphene, 3.0), (square, 4.0), (cubic, 6.0)])
def test_the_interval_contains_the_spectrum(build, half_width):
    """The property that matters. Too wide costs resolution; too narrow is fatal."""
    lo, hi = kite.estimate_bounds(build())
    assert lo < -half_width and hi > half_width


def test_the_margin_is_generous_but_not_wasteful():
    """10% headroom by PM decision (2026-09-16).

    The two failure modes are not symmetric: under-padding puts eigenvalues
    outside [-1, 1], where the run still exits 0 and writes meaningless numbers,
    while over-padding only costs resolution, linearly. A small amount of
    Anderson disorder the caller forgot to declare is enough to cross a 1%
    margin, so 1% was too tight.
    """
    lo, hi = kite.estimate_bounds(graphene())
    padding = (hi - lo) / 6.0
    assert padding == pytest.approx(1.10, rel=1e-6)
    assert padding < 1.5, "but it must not throw away resolution either"


def test_uniform_disorder_uses_the_standard_deviation_correctly():
    """sigma is a standard deviation; the box half-width is sqrt(3) times larger."""
    lat = graphene()
    disorder = kite.Disorder(lat)
    disorder.add_disorder("A", "Uniform", 0.0, 0.2)
    lo, hi = kite.estimate_bounds(lat, disorder=disorder)
    assert hi == pytest.approx((3.0 + np.sqrt(3) * 0.2) * 1.10, rel=1e-6)


def test_gaussian_disorder_uses_a_tail_cut():
    lat = graphene()
    disorder = kite.Disorder(lat)
    disorder.add_disorder("A", "Gaussian", 0.0, 0.1)
    lo, hi = kite.estimate_bounds(lat, disorder=disorder)
    assert hi == pytest.approx((3.0 + 8.0 * 0.1) * 1.10, rel=1e-6)


def test_a_magnetic_field_falls_back_to_gershgorin():
    """Hofstadter broadening exceeds the zero-field band edges, so the BZ scan
    is not a bound once a field is applied."""
    lat = graphene()
    plain = kite.estimate_bounds(lat)[1]
    with_field = kite.estimate_bounds(lat, magnetic_field=True)[1]
    assert with_field >= plain


def test_it_beats_the_automatic_estimate_for_graphene():
    """The automatic path uses max_onsite + 6*max_hopping, which for a
    three-coordinated lattice is twice the true half-bandwidth -- and the
    resolution sigma = pi*energy_scale/M degrades in direct proportion."""
    lo, hi = kite.estimate_bounds(graphene(t=2.8))
    automatic = 6 * 2.8
    assert hi < 0.6 * automatic
    assert hi > 3 * 2.8, "but it must still contain the spectrum"


# ---------------------------------------------------------------------------
# The automatic path inside config_system (T-023)
# ---------------------------------------------------------------------------

def _export(lattice, tmp_path, **kwargs):
    import h5py
    dim = len(lattice.vectors)                 # 2D and 3D lattices need matching lengths
    configuration = kite.Configuration(
        divisions=[1] * dim, length=[16] * dim, boundaries=["periodic"] * dim,
        is_complex=bool(kwargs.get("modification")), precision=1)   # no spectrum_range
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=100, num_moments=32, num_random=1)
    path = str(tmp_path / "auto.h5")
    kite.config_system(lattice, configuration, calculation, filename=path, **kwargs)
    with h5py.File(path, "r") as f:
        return float(f["EnergyScale"][()])


def test_automatic_scaling_uses_the_estimator(tmp_path):
    """Omitting spectrum_range must give the real half-bandwidth plus margin,
    not the coordination-six rule.

    For graphene the old expression returned 6*|t| against a true 3*|t| -- twice
    too wide -- and since sigma = pi*energy_scale/M that halved the resolution
    every user got by default.
    """
    scale = _export(graphene(t=2.8), tmp_path)
    assert scale == pytest.approx(3 * 2.8 * 1.10, rel=1e-6)
    assert scale < 6 * 2.8, "must be tighter than the old coordination-six rule"


def test_automatic_scaling_still_contains_the_spectrum(tmp_path):
    for build, half_width in ((graphene, 3.0), (square, 4.0), (cubic, 6.0)):
        assert _export(build(), tmp_path) > half_width


def test_automatic_scaling_accounts_for_disorder(tmp_path):
    """The old rule ignored disorder entirely."""
    lat = graphene()
    disorder = kite.Disorder(lat)
    disorder.add_disorder("A", "Uniform", 0.0, 0.3)
    clean = _export(graphene(), tmp_path)
    dirty = _export(lat, tmp_path, disorder=disorder)
    assert dirty > clean


def test_a_magnetic_field_does_not_use_the_band_edges(tmp_path):
    """Hofstadter broadening exceeds the zero-field edges, so the BZ scan is
    not a bound once a field is applied."""
    scale = _export(graphene(), tmp_path, modification=kite.Modification(flux=1 / 16))
    assert scale >= 3.0, "must still contain the zero-field spectrum at least"
