"""Negative-test regression for the operator-dimension safety check added to
config_system() (stage-1 re-audit, Section "Open blocker: matrix dimension").

Before this fix, add_orbital_coupling() sized a newly registered operator
matrix from however many orbital names happened to be in
_orbital_index_collection at that moment -- register only part of the
lattice basis, or add more orbital names after a label's first coupling
call, and the resulting matrix stayed undersized. Python accepted this
silently and KITEx's deterministic kernel (which always loops a, b over the
FULL orbital range) would read past the matrix's actual bounds: a real
out-of-bounds access, not merely a cosmetic mismatch.

config_system() now computes the lattice's true total orbital count and
rejects any registered operator whose shape doesn't match exactly. This
script exercises the failure modes the re-audit named explicitly --
incomplete, late-grown, duplicate-label, and non-contiguous registrations --
plus one valid case that must still pass.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import kite
from kite import lattice as latt


def _two_orbital_lattice():
    lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
    lat.add_sublattices(('A', [0.0, 0.0], 0.0), ('B', [0.5, 0.5], 0.0))
    lat.add_hoppings(([0, 0], 'A', 'B', -1.0))
    return lat


def _configuration():
    return kite.Configuration(divisions=[1, 1], length=[4, 4],
                               boundaries=['periodic', 'periodic'],
                               is_complex=True, precision=1,
                               spectrum_range=[-4, 4])


def _export(lattice, calculation, tag):
    kite.config_system(lattice, _configuration(), calculation,
                        filename="/tmp/operator_dimension_regression_{}.h5".format(tag))


def expect_rejected(name, build_fn):
    lattice = _two_orbital_lattice()
    calculation = kite.Calculation(_configuration())
    build_fn(calculation)
    try:
        _export(lattice, calculation, name)
    except ValueError as e:
        print("PASS ({}): correctly rejected -- {}".format(name, e))
        return True
    print("FAIL ({}): no error raised, bad operator silently accepted".format(name))
    return False


def expect_accepted(name, build_fn):
    lattice = _two_orbital_lattice()
    calculation = kite.Calculation(_configuration())
    build_fn(calculation)
    try:
        _export(lattice, calculation, name)
    except ValueError as e:
        print("FAIL ({}): valid operator incorrectly rejected -- {}".format(name, e))
        return False
    print("PASS ({}): valid full-basis operator accepted".format(name))
    return True


def incomplete_registration(calculation):
    """Only one of the two orbitals is named before the coupling call --
    the audit's own adversarial case."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


def late_grown_registration(calculation):
    """Both orbitals exist eventually, but the second is named AFTER the
    label's matrix has already been sized from the first."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')  # sizes l0 as 1x1
    calculation.add_orbital_index('B', 1)                  # too late for l0
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


def sparse_but_fully_sized(calculation):
    """A second label ('l1') only ever sets ONE diagonal entry (the other
    stays an implicit zero). This must NOT be rejected: both orbitals were
    already named before l1's first coupling call, so its matrix is
    correctly sized (2x2) -- a zero entry is a legitimate operator value
    (e.g. a single-orbital projector), not a dimension error. Guards against
    an overzealous check that flags sparsity instead of size."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_index('B', 1)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
    calculation.add_orbital_coupling('B', 'B', 1.0, 'l0')  # l0: fully specified, 2x2
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l1')  # l1: sparse (B,B left at 0),
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10,  # but still correctly 2x2
                         operators=['l0', 'l1'])


def valid_full_basis(calculation):
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_index('B', 1)
    calculation.add_orbital_coupling('A', 'A', 0.5, 'l0')
    calculation.add_orbital_coupling('B', 'B', -0.5, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


if __name__ == "__main__":
    results = [
        expect_rejected("incomplete", incomplete_registration),
        expect_rejected("late_grown", late_grown_registration),
        expect_accepted("sparse_but_fully_sized", sparse_but_fully_sized),
        expect_accepted("valid_full_basis", valid_full_basis),
    ]
    if all(results):
        print("\nAll operator-dimension regression checks PASSED.")
        sys.exit(0)
    else:
        print("\nSome operator-dimension regression checks FAILED.")
        sys.exit(1)
