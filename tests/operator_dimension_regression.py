"""Negative-test regression for config_system()'s operator/orbital-index
safety checks (stage-1 re-audit and third-pass audit).

Three distinct defects, three distinct checks:

1. Undersized operator matrix (re-audit, Section "Open blocker: matrix
   dimension"): add_orbital_coupling() sizes a newly registered operator
   matrix from however many orbital names happened to be in
   _orbital_index_collection at that moment -- register only part of the
   lattice basis, or add more orbital names after a label's first coupling
   call, and the resulting matrix stayed undersized. KITEx's deterministic
   kernel (which always loops a, b over the FULL orbital range) would read
   past the matrix's actual bounds. Fixed by checking every operator's shape
   against the lattice's true total orbital count at export time.

2. Duplicate/non-contiguous orbital indices (third-pass audit, Section "Open
   API footgun: orbital-index aliases"): add_orbital_index() accepted any
   integer index with no uniqueness or range check. Two names mapped to the
   same index still produce a CORRECTLY SIZED matrix (so check #1 above
   can't catch it) but couplings targeting either name silently overwrite
   the same entry. Fixed by requiring the index map to be a bijection onto
   0..N_orb-1.

3. Post-registration Hermiticity mutation (third-pass audit, Section "Open
   blocker: validation is not final"): ldos()/ldos_map() validate
   Hermiticity when called, but a later add_orbital_coupling() call on the
   same label can mutate it into a non-Hermitian matrix before export, and
   the original check has no way to see that. Fixed by re-validating the
   final matrix, at export time, for every label actually referenced by an
   ldos()/ldos_map() request.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
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
    with tempfile.TemporaryDirectory(prefix="kite_operator_dimension_regression_") as tmpdir:
        kite.config_system(lattice, _configuration(), calculation,
                            filename=os.path.join(tmpdir, "{}.h5".format(tag)))


def expect_rejected(name, build_fn):
    """Some malformed registrations are caught by config_system()'s
    export-time checks (ValueError); others -- any out-of-range index whose
    matrix position exceeds the current (name-count-sized) matrix -- are
    already caught earlier, by plain numpy bounds checking inside
    add_orbital_coupling() itself (IndexError). Both are legitimate
    "rejected", since either prevents the malformed operator from silently
    reaching KITEx."""
    lattice = _two_orbital_lattice()
    calculation = kite.Calculation(_configuration())
    try:
        build_fn(calculation)
        _export(lattice, calculation, name)
    except (ValueError, IndexError) as e:
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


def duplicate_index_registration(calculation):
    """Two orbital names aliased to the same index -- the third-pass audit's
    own reproduction. Matrix ends up correctly SIZED (2x2, since two names
    were registered), so check #1 above can't catch this; couplings to A and
    B silently overwrite the same (0,0) entry instead of writing two
    different diagonal entries."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_index('B', 0)
    calculation.add_orbital_coupling('A', 'A', -1.0, 'l0')
    calculation.add_orbital_coupling('B', 'B', 1.0, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


def non_contiguous_registration(calculation):
    """Indices {0, 2} for a 2-orbital lattice -- valid count, valid
    uniqueness, but skips index 1 and reaches out of range at index 2."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_index('B', 2)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
    calculation.add_orbital_coupling('B', 'B', -1.0, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


def negative_index_registration(calculation):
    """A negative index is unique and in a sense "contiguous" with 0, but
    not a valid position in a 0-based orbital basis."""
    calculation.add_orbital_index('A', -1)
    calculation.add_orbital_index('B', 1)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
    calculation.add_orbital_coupling('B', 'B', -1.0, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])


def post_registration_hermiticity_mutation(calculation):
    """l0 is Hermitian (diagonal, real) when ldos_map() validates it, but a
    later add_orbital_coupling() call adds an off-diagonal entry with no
    Hermitian-conjugate partner, mutating it into a non-Hermitian matrix
    before export -- the third-pass audit's own reproduction."""
    calculation.add_orbital_index('A', 0)
    calculation.add_orbital_index('B', 1)
    calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
    calculation.ldos_map(energy_=0.0, sigma_=0.1, vectors_=10, operators=['l0'])
    calculation.add_orbital_coupling('A', 'B', 1j, 'l0')  # no conjugate partner


if __name__ == "__main__":
    results = [
        expect_rejected("incomplete", incomplete_registration),
        expect_rejected("late_grown", late_grown_registration),
        expect_accepted("sparse_but_fully_sized", sparse_but_fully_sized),
        expect_accepted("valid_full_basis", valid_full_basis),
        expect_rejected("duplicate_index", duplicate_index_registration),
        expect_rejected("non_contiguous", non_contiguous_registration),
        expect_rejected("negative_index", negative_index_registration),
        expect_rejected("post_registration_hermiticity_mutation",
                        post_registration_hermiticity_mutation),
    ]
    if all(results):
        print("\nAll operator-dimension regression checks PASSED.")
        sys.exit(0)
    else:
        print("\nSome operator-dimension regression checks FAILED.")
        sys.exit(1)
