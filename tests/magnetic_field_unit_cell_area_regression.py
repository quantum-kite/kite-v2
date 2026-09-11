"""Regression for the 2D unit-cell area used by the automated magnetic field.

Defect: config_system() computed the unit-cell area with

    np.linalg.norm(np.cross(vectors[0, :], vectors[1, :]))

NumPy 2.0 removed the two-dimensional form of np.cross, so for any 2D lattice
this raised

    ValueError: Both input arrays must be (arrays of) 3-dimensional vectors

The call sits in the `if modification.magnetic_field or modification.flux:`
branch, before either sub-branch, so BOTH entry points died -- every 2D
magnetic-field calculation, on any NumPy >= 2.0. 3D was unaffected, because
three-element vectors still have a cross product.

The fix embeds both lattice vectors in 3D before crossing. For a 2D cell only
the z-component of the product survives, |a1_x a2_y - a1_y a2_x|, which is the
area; for a 3D cell the padding is a copy and the arithmetic is unchanged.

Two things are asserted here, and the second matters as much as the first:

1. The 2D area is correct -- graphene with a = 1 has |a1 x a2| = sqrt(3)/2.
2. The 3D result is bit-identical to the pre-fix expression. A fix that
   silently perturbed 3D would be a regression in a working capability, and
   floating-point equality is the right bar because the padded operands are
   the same numbers in the same order.
"""

import numpy as np
import pytest


def _cell_area(vectors):
    """The expression now used by config_system(), isolated for test."""
    a1 = np.zeros(3)
    a2 = np.zeros(3)
    a1[:vectors.shape[1]] = vectors[0, :]
    a2[:vectors.shape[1]] = vectors[1, :]
    return np.linalg.norm(np.cross(a1, a2))


def _cell_area_prefix(vectors):
    """The pre-fix expression. Valid only for 3-vectors; kept as the oracle."""
    return np.linalg.norm(np.cross(vectors[0, :], vectors[1, :]))


def test_two_dimensional_area_is_computed_at_all():
    """The defect itself: this raised ValueError on NumPy >= 2.0."""
    graphene = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    assert _cell_area(graphene) == pytest.approx(np.sqrt(3) / 2)


def test_two_dimensional_area_matches_the_determinant():
    """|a1 x a2| in 2D is |det[a1; a2]|, for a few unrelated cells."""
    for cell in (
        np.array([[1.0, 0.0], [0.0, 1.0]]),           # square
        np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]]),  # hexagonal
        np.array([[2.0, 0.3], [-0.4, 1.7]]),          # oblique
    ):
        assert _cell_area(cell) == pytest.approx(abs(np.linalg.det(cell)))


def test_three_dimensional_result_is_bit_identical():
    """A working capability must not move, not even in the last bit."""
    for cell in (
        np.eye(3),
        np.array([[1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0], [0.0, 0.0, 2.0]]),
        np.array([[1.3, -0.2, 0.1], [0.4, 2.1, -0.7], [0.0, 0.0, 3.2]]),
    ):
        assert _cell_area(cell) == _cell_area_prefix(cell)


def test_area_is_positive_for_a_left_handed_cell():
    """Ordering of the basis vectors must not flip the sign of an area."""
    right = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    left = np.array([[0.5, np.sqrt(3) / 2], [1.0, 0.0]])
    assert _cell_area(left) > 0
    assert _cell_area(left) == pytest.approx(_cell_area(right))
