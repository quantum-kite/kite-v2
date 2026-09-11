"""Negative-test regression for the first-neighbour hopping-range limit.

KITE packs a hopping's destination into a single integer:

    orbital * 3**D + sum_d (n_d + 1) * 3**d

so the cell field is exactly D base-3 digits -- -1, 0, +1, and nothing else.
A component of 2 does not reach a further cell. It carries into the digit
above: for any direction but the last that is the next CELL digit (wrong cell,
valid orbital); for the last direction it is the ORBITAL field, and the C++
decoder then indexes past the orbital count.

Either way the file is well formed, the run is clean, and the spectrum belongs
to a different Hamiltonian. That is the worst failure mode available -- a wrong
answer with no signal -- which is why the limit is enforced at export.

The constraint is a property of the file format, confirmed with the KITE
authors as intentional. Raising NGHOSTS does not relax it: NGHOSTS is the halo
depth, so the neighbouring cell's data is in memory, but the address cannot
name it.

The identical check already existed on the structural-disorder path
(`add_local_bond_disorder`) and was missing on the Hamiltonian path, which is
the one that matters more.
"""

import numpy as np
import pytest

import kite
from kite import lattice as latt


def _square(relative_index, t_extra=-0.3):
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(("A", [0, 0], 0.0))
    lat.add_hoppings(([1, 0], "A", "A", -1.0),
                     (relative_index, "A", "A", t_extra))
    return lat


def _export(lattice, tmp_path):
    configuration = kite.Configuration(
        divisions=[1, 1], length=[16, 16],
        boundaries=["periodic", "periodic"],
        is_complex=False, precision=1, spectrum_range=[-6, 6])
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=100, num_moments=32, num_random=1, num_disorder=1)
    return kite.config_system(lattice, configuration, calculation,
                              filename=str(tmp_path / "probe.h5"))


@pytest.mark.parametrize("relative_index", [[0, 1], [1, 1], [-1, 1], [-1, -1], [1, 0]])
def test_first_neighbour_cells_are_accepted(relative_index, tmp_path):
    """Every corner of the 3x3 block is representable and must export."""
    _export(_square(relative_index), tmp_path)


@pytest.mark.parametrize("relative_index", [[2, 0], [0, 2], [-2, 0], [1, -2], [2, 2], [3, 1]])
def test_beyond_first_neighbour_is_refused(relative_index, tmp_path):
    """Anything outside -1..+1 must be refused, not silently re-encoded."""
    with pytest.raises(SystemExit) as excinfo:
        _export(_square(relative_index), tmp_path)
    message = str(excinfo.value)
    assert "out of range" in message
    assert str(relative_index) in message, "the message must name the offending hopping"


def test_the_shipped_paper_model_is_refused(tmp_path):
    """examples/paper/Section_6_benchmarks/Input/dos_gr_nnn.py:62 declares

        [(-2, +1), 'A', 'B', t_nnn]

    on a honeycomb lattice -- a real, shipped model that produces a wrong
    Hamiltonian today. The script itself still needs pybinding, so the bond is
    reproduced natively here.
    """
    a_cc = 1 / np.sqrt(3)
    lat = latt.Lattice(a1=[1, 0], a2=[0.5, np.sqrt(3) / 2])
    lat.add_sublattices(("A", [0, -a_cc / 2], 0.0), ("B", [0, +a_cc / 2], 0.0))
    lat.add_hoppings(([0, 0], "A", "B", -2.8),
                     ([-1, 0], "A", "B", -2.8),
                     ([-1, 1], "A", "B", -2.8),
                     ([-2, +1], "A", "B", -0.1))      # the violating next-nearest bond

    with pytest.raises(SystemExit) as excinfo:
        _export(lat, tmp_path)
    assert "[-2, 1]" in str(excinfo.value)


def test_the_message_does_not_offer_nghosts_as_an_escape(tmp_path):
    """docs/documentation/tb_model.md offered NGHOSTS as a way round the limit.
    It is not one, and the error must not repeat the claim."""
    with pytest.raises(SystemExit) as excinfo:
        _export(_square([2, 0]), tmp_path)
    message = str(excinfo.value)
    assert "NGHOSTS does not help" in message
    assert "larger unit cell" in message, "the message should say what to do instead"
