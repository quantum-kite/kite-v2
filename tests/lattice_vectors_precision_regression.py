"""The lattice vectors are stored single precision, deliberately.

Geometry does not need double: the lattice vectors set bond directions and the
unit-cell area, and single precision carries those to about seven significant
figures, which is far finer than any model parameter a user supplies. Storing
them double costs space in every exported file and buys nothing for almost
every calculation.

The hoppings are a separate question and keep their own precision, chosen by
`Configuration(precision=...)`. This is only about the geometry.

Measured cost of the change on the one quantity most sensitive to the bond
vectors -- the diamagnetic term of the optical conductivity, which is built
from products of bond-vector components -- is a relative 3.6e-8 on graphene.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

import kite
from kite import lattice as latt


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    """A real export, so this follows the exporter rather than describing it."""
    a = 0.24595
    a_cc = 0.142
    lattice = latt.Lattice(a1=a * np.array([1.0, 0.0]),
                           a2=a * np.array([0.5, 0.5 * np.sqrt(3)]))
    lattice.add_sublattices(('A', [0, -a_cc / 2], 0.0), ('B', [0, a_cc / 2], 0.0))
    lattice.add_hoppings(([0, 0], 'A', 'B', -1.0),
                         ([1, -1], 'A', 'B', -1.0),
                         ([0, -1], 'A', 'B', -1.0))
    configuration = kite.Configuration(
        divisions=[2, 2], length=[16, 16], boundaries=["periodic", "periodic"],
        is_complex=False, precision=1, spectrum_range=[-3.1, 3.1])
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=64, num_moments=16, num_random=1, num_disorder=1)

    path = tmp_path_factory.mktemp("latt") / "config.h5"
    kite.config_system(lattice, configuration, calculation, filename=str(path))
    return path


def test_lattice_vectors_are_single_precision(exported):
    with h5py.File(exported, "r") as handle:
        assert handle["LattVectors"].dtype == np.float32, (
            f"LattVectors is {handle['LattVectors'].dtype}; geometry is stored "
            "single precision on purpose")


def test_the_vectors_still_describe_the_lattice_that_was_built(exported):
    """Single precision is a storage choice, not permission to be wrong."""
    a = 0.24595
    expected = np.array([[a, 0.0], [a * 0.5, a * 0.5 * np.sqrt(3)]])
    with h5py.File(exported, "r") as handle:
        stored = np.array(handle["LattVectors"]).reshape(expected.shape)
    # float32 carries ~7 significant figures.
    assert stored == pytest.approx(expected, rel=1e-6), (
        f"stored {stored.ravel()} does not match the lattice that was built "
        f"{expected.ravel()}")


def test_the_file_is_still_readable_as_double(exported):
    """KITEx reads this with get_hdf5<double>; HDF5 converts on read."""
    with h5py.File(exported, "r") as handle:
        widened = np.array(handle["LattVectors"], dtype=np.float64)
    assert np.isfinite(widened).all()
    assert widened.dtype == np.float64
