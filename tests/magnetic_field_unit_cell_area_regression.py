"""Regression for the 2D unit-cell area used by the automated magnetic field.

Defect: `config_system()` computed the unit-cell area with

    np.linalg.norm(np.cross(vectors[0, :], vectors[1, :]))

NumPy 2.0 removed the two-dimensional form of `np.cross`, so for any 2D lattice
this raised

    ValueError: Both input arrays must be (arrays of) 3-dimensional vectors

The call sits in the `if modification.magnetic_field or modification.flux:`
branch, before either sub-branch, so BOTH entry points died -- every 2D
magnetic-field calculation, on any NumPy >= 2.0. 3D was unaffected, because
three-element vectors still have a cross product. The fix embeds both lattice
vectors in 3D before crossing: for a 2D cell only the z-component survives,
|a1_x a2_y - a1_y a2_x|, which is the area; for 3D the padding is a copy.

**This file drives `config_system` and reads what it exported.** It used to
hold a hand copy of the production expression and assert things about the copy,
which is no test of the shipped code at all: the copy had drifted already --
it omitted the `* 1e-18` that converts nm^2 to m^2 -- and deleting the fix from
`python/kite/__init__.py` left every case passing. The observable here is
`/Hamiltonian/MagneticFieldMul`, which is
`round(magnetic_field / magnetic_field_min)` with

    magnetic_field_min = phi0 / (length[1] * unit_cell_area)

so a unit-cell area that is wrong by any factor moves the exported integer.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

import kite
from kite import lattice as latt

#: As `config_system` defines them. Repeated here because the expected field is
#: the point of comparison -- these are physical constants, not an expression
#: lifted out of the code under test.
HBAR_EV_S = 6.58211899e-16
PHI0 = 2 * np.pi * HBAR_EV_S
NM2_TO_M2 = 1e-18


def _expected_minimum_field(vectors, length_y):
    """phi0 / (L_y * A), with A the cell area in m^2, from the geometry alone."""
    a1, a2 = np.zeros(3), np.zeros(3)
    a1[:vectors.shape[1]] = vectors[0]
    a2[:vectors.shape[1]] = vectors[1]
    area = abs(np.linalg.det(np.array([a1[:2], a2[:2]]))) if vectors.shape[1] == 2 \
        else np.linalg.norm(np.cross(a1, a2))
    return PHI0 / (length_y * area * NM2_TO_M2)


def _square_2d(a=1.0):
    lat = latt.Lattice(a1=[a, 0], a2=[0, a])
    lat.add_sublattices(("A", [0, 0], 0.0))
    lat.add_hoppings(([1, 0], "A", "A", -1.0), ([0, 1], "A", "A", -1.0))
    return lat


def _hexagonal_2d(a=1.0):
    lat = latt.Lattice(a1=[a, 0], a2=[0.5 * a, np.sqrt(3) / 2 * a])
    lat.add_sublattices(("A", [0, 0], 0.0), ("B", [0, a / np.sqrt(3)], 0.0))
    lat.add_hoppings(([0, 0], "A", "B", -1.0), ([1, -1], "A", "B", -1.0),
                     ([0, -1], "A", "B", -1.0))
    return lat


def _cubic_3d(a=1.0):
    lat = latt.Lattice(a1=[a, 0, 0], a2=[0, a, 0], a3=[0, 0, a])
    lat.add_sublattices(("A", [0, 0, 0], 0.0))
    lat.add_hoppings(([1, 0, 0], "A", "A", -1.0), ([0, 1, 0], "A", "A", -1.0),
                     ([0, 0, 1], "A", "A", -1.0))
    return lat


def _export(lattice, length, field, tmp_path, divisions=None):
    dimensions = len(length)
    configuration = kite.Configuration(
        divisions=divisions or [1] * dimensions, length=length,
        boundaries=["periodic"] * dimensions, is_complex=True, precision=1,
        spectrum_range=[-10.0, 10.0])
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=64, num_moments=32, num_random=1, num_disorder=1)
    path = tmp_path / "field.h5"
    kite.config_system(lattice, configuration, calculation, filename=str(path),
                       modification=kite.Modification(magnetic_field=field))
    with h5py.File(path, "r") as handle:
        return float(np.array(handle["/Hamiltonian/MagneticFieldMul"]))


@pytest.mark.parametrize("build,vectors", [
    (_square_2d, np.array([[1.0, 0.0], [0.0, 1.0]])),
    (_hexagonal_2d, np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])),
])
def test_a_two_dimensional_lattice_exports_at_all(build, vectors, tmp_path):
    """The defect itself: this raised ValueError on NumPy >= 2.0."""
    minimum = _expected_minimum_field(vectors, 64)
    assert _export(build(), [64, 64], minimum, tmp_path) == pytest.approx(1.0)


@pytest.mark.parametrize("multiple", [1, 3, 7])
def test_the_exported_multiple_follows_the_cell_area(multiple, tmp_path):
    """Asking for N minimum fields must export N, for the real cell area."""
    vectors = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    minimum = _expected_minimum_field(vectors, 64)
    exported = _export(_hexagonal_2d(), [64, 64], multiple * minimum, tmp_path)
    assert exported == pytest.approx(float(multiple)), (
        f"asked for {multiple} x the minimum commensurate field and the export "
        f"records {exported}; config_system's unit-cell area disagrees with "
        "|a1 x a2| for this lattice")


def test_the_area_carries_the_lattice_constant(tmp_path):
    """Halving `a` quarters the cell, so the minimum field quadruples.

    The hand-copied expression this file used to test omitted the nm^2 -> m^2
    conversion, which no comparison against another copy could notice. A ratio
    taken from two real exports does.
    """
    big = _expected_minimum_field(np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]]), 64)
    small = _expected_minimum_field(
        np.array([[0.5, 0.0], [0.25, np.sqrt(3) / 4]]), 64)
    assert small / big == pytest.approx(4.0)
    assert _export(_hexagonal_2d(a=0.5), [64, 64], small, tmp_path) == pytest.approx(1.0)


def test_three_dimensional_still_works(tmp_path):
    """A working capability must not move: 3D never had the defect."""
    vectors = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    minimum = _expected_minimum_field(vectors, 32)
    assert _export(_cubic_3d(), [32, 32, 32], 5 * minimum, tmp_path) == pytest.approx(5.0)


def test_the_basis_order_does_not_flip_the_area(tmp_path):
    """|a1 x a2| is an area; swapping the vectors must not change the field."""
    right = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    left = np.array([[0.5, np.sqrt(3) / 2], [1.0, 0.0]])
    assert _expected_minimum_field(left, 64) == pytest.approx(
        _expected_minimum_field(right, 64))

    lattice = latt.Lattice(a1=[0.5, np.sqrt(3) / 2], a2=[1.0, 0.0])
    lattice.add_sublattices(("A", [0, 0], 0.0))
    lattice.add_hoppings(([1, 0], "A", "A", -1.0), ([0, 1], "A", "A", -1.0))
    minimum = _expected_minimum_field(left, 64)
    assert _export(lattice, [64, 64], minimum, tmp_path) == pytest.approx(1.0)
