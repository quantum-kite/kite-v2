"""KITE-tools must decode the spectral_map k index, and decode it correctly.

`/Calculation/spectral_map/Map` is a flat array over the reciprocal grid, and
`store_spectral` re-indexes it by `(coord + L/2) % L` on every axis before
writing. So the reader has two things to get right -- the flat ordering and the
shift -- and getting either wrong produces a map that still looks like a map,
with the weight on the wrong k-points.

The reference is a dispersion we know in closed form. On a lattice of
nearly-decoupled chains (`t` along one lattice vector, 1/20 of it along the
other) the band is essentially `E(k) = -2t cos(k.a1)`, so at `E = -2t cos(pi/4)`
the weight sits on `k.a1 = +-pi/4` -- that is, at fractional coordinate
`+-1/8` -- for every value of the other component.

`pi/4` is chosen, not incidental. `k = +-pi/2` is a fixed point of the
`(s + L/2) % L` re-indexing, so a contour sitting there cannot tell a shifted
layout from an unshifted one; `pi/4` can, and an unshifted reader puts the
weight at `+-3/8` instead.

The lattice is triangular rather than square so that the Cartesian conversion
is not the identity. Reading `k.a1 = +-pi/4` off the Cartesian columns requires
the reciprocal vectors to be `2*pi*vectors^{-T}`; with `vectors^{-1}` the two
reciprocal vectors are exchanged and `kx` is no longer constant along the
contour.
"""

import re
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

LX, LY = 16, 32          # multiples of TILE=8, and different, so no ordering ties
ENERGY = -2 * np.cos(np.pi / 4)
SIGMA = 2e-2

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

strong_axis = int(sys.argv[2])          # 0: hopping along a1, 1: along a2
t = [0.05, 0.05]
t[strong_axis] = 1.0

# triangular, so the reciprocal vectors are not parallel to the Cartesian axes
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.5, np.sqrt(3) / 2])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', -t[0]), ([0, 1], 'A', 'A', -t[1]))

cfg = kite.Configuration(divisions=[1, 1], length=[16, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-3, 3])
calc = kite.Calculation(cfg)
calc.spectral_map(energy_=-2 * np.cos(np.pi / 4), sigma_=2e-2, vectors_=16)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _data(path):
    rows = [line.split() for line in path.read_text().splitlines()
            if not line.startswith("#")]
    assert rows, f"{path.name} has no data rows"
    return np.array([[float(c) for c in r] for r in rows])


def _contour(data):
    """The rows carrying the spectral weight: columns f1, f2, kx, ky."""
    value = data[:, 5]
    return data[value > 0.1 * value.max()]


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    out = {}
    for axis in (0, 1):
        d = tmp_path_factory.mktemp(f"spectralmap{axis}")
        subprocess.run([sys.executable, "-c", _BUILDER, "m.h5", str(axis)],
                       cwd=d, check=True, capture_output=True, text=True)
        subprocess.run([str(KITEX), "m.h5"], cwd=d, check=True, capture_output=True)
        subprocess.run([str(KITE_TOOLS), "m.h5", "--SpectralMap", "-X", "-N", "m.dat"],
                       cwd=d, check=True, capture_output=True)
        out[axis] = _data(d / "m.dat")
    return out


def test_every_k_point_is_written_once(runs):
    data = runs[0]
    assert len(data) == LX * LY, (
        f"expected one row per k-point per orbital, got {len(data)}")
    grid = {(round(f1, 6), round(f2, 6)) for f1, f2 in data[:, :2]}
    assert len(grid) == LX * LY, "the decoded k-points are not distinct"
    assert min(f for f, _ in grid) == pytest.approx(-0.5)
    assert max(f for f, _ in grid) == pytest.approx(0.5 - 1.0 / LX)


def test_the_fftshift_is_undone(runs):
    """The contour must land on k.a1 = +-pi/4, i.e. fractional +-1/8."""
    for axis in (0, 1):
        contour = _contour(runs[axis])
        along = sorted(set(np.round(contour[:, axis], 6)))
        assert along == [-0.125, 0.125], (
            f"with the strong hopping along a{axis + 1} the weight decoded to "
            f"fractional {along}, not [-0.125, 0.125]. store_spectral shifts "
            "every axis by (s + L/2) % L; leaving that shift in place puts the "
            "weight at +-0.375")


def test_the_two_axes_are_not_interchangeable(runs):
    """Moving the strong hopping must move the contour to the other axis."""
    for axis, other in ((0, 1), (1, 0)):
        contour = _contour(runs[axis])
        spread = contour[:, other].max() - contour[:, other].min()
        assert spread > 0.5, (
            f"with the strong hopping along a{axis + 1} the contour is confined "
            f"in the a{other + 1} direction too, so this test could not tell the "
            "two axes apart")


def test_the_cartesian_components_use_the_reciprocal_basis(runs):
    """k.a1 = 2*pi*f1, and a1 is the x axis, so kx is pi/4 on the contour."""
    contour = _contour(runs[0])
    kx = sorted(set(np.round(contour[:, 2], 6)))
    assert kx == [pytest.approx(-np.pi / 4), pytest.approx(np.pi / 4)], (
        f"kx on the contour is {kx}, not +-pi/4. The reciprocal basis is "
        "2*pi*vectors^{-T}; with vectors^{-1} the two reciprocal vectors are "
        "exchanged")
    ky_spread = contour[:, 3].max() - contour[:, 3].min()
    assert ky_spread > 1.0, (
        "ky does not vary along the contour, so a conversion that ignored the "
        "second reciprocal vector would pass this test")


def test_the_energy_is_written_back_in_eV(runs, tmp_path_factory):
    """Energy and width are stored rescaled; the header must undo that."""
    path = next(iter(tmp_path_factory.getbasetemp().glob("spectralmap0*/m.dat")))
    header = path.read_text().splitlines()[0]
    match = re.search(r"at E = ([-\d.eE+]+) eV, sigma = ([-\d.eE+]+) eV", header)
    assert match, f"the header no longer states E and sigma in eV: {header!r}"
    assert float(match.group(1)) == pytest.approx(ENERGY, rel=1e-5)
    assert float(match.group(2)) == pytest.approx(SIGMA, rel=1e-5)


def test_the_absolute_normalisation_is_the_spectral_function_in_1_per_eV(runs):
    """The header claims 1/eV, so check the values against the closed form.

    One orbital and no disorder, so the spectral function is exactly
    `g_sigma(E - eps_k)` with `eps_k = -2 t1 cos(2 pi f1) - 2 t2 cos(2 pi f2)`
    and `g` a unit-area Gaussian in eV. This is the one claim in the header a
    user can check, and it checks the decode at the same time: a wrong k label
    would compare each value against the wrong `eps_k`.
    """
    data = runs[0]
    f1, f2, value = data[:, 0], data[:, 1], data[:, 5]
    eps = -2 * 1.0 * np.cos(2 * np.pi * f1) - 2 * 0.05 * np.cos(2 * np.pi * f2)
    exact = np.exp(-(ENERGY - eps) ** 2 / (2 * SIGMA ** 2)) / (SIGMA * np.sqrt(2 * np.pi))

    assert np.corrcoef(value, exact)[0, 1] > 1 - 1e-6, (
        "the map does not follow g_sigma(E - eps_k) over the grid")
    assert value.sum() == pytest.approx(exact.sum(), rel=1e-3), (
        "the map is the right shape but the wrong size: the gaussian branch is "
        "normalised so that integrating over E gives 1 for each (k, orbital), "
        "which fixes the overall scale and makes 1/eV the honest unit")


def test_the_units_are_stated(runs, tmp_path_factory):
    path = next(iter(tmp_path_factory.getbasetemp().glob("spectralmap0*/m.dat")))
    text = path.read_text()
    assert "1/eV" in text, "the gaussian branch no longer states its units"
    assert "fractional reciprocal-lattice" in text, (
        "the header no longer says which of the two k columns is which")
    assert "k.a_d = 2*pi*f_d" in text, (
        "the header no longer says where the 2*pi lives; gaussian_wave_packet "
        "takes k in the opposite convention, so this is worth stating")
    assert "fixed-twist" in text, (
        "the header no longer warns that random boundary twists smear the k grid")
