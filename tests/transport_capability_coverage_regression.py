"""Each transport capability names the test that protects it, and that test passes.

RF-007 to RF-009 are the conductivities -- Kubo-Bastin DC, single-shot DC and the
AC optical response -- and the legacy suite covers all three, on two different
lattices each so that a 2D and a 3D case both run.

RF-011 and RF-012 are the vertex-operator language, and they are covered by the
modern tests rather than the legacy suite, which predates it. RF-012 is the one
that is *not* fully discharged: determining sign and normalisation from the
user's vertex definition rather than from a runtime heuristic is gated on a
demonstrated failing case (DEC-013), so until then what is protected is the
documented heuristic itself -- the `NumVelocities mod 4` rule and the units the
readers print. That is what the tests named below pin, and naming them here is
not a claim that the requirement is satisfied.
"""

import pytest

from capability_coverage import (BINARIES_BUILT, TESTS,
                                 assert_the_directory_passes, named_tests_exist)

#: capability -> (legacy directories, modern test files)
COVERAGE = {
    "RF-007": (["test009_KITEx_sq_dc", "test016_KITEx_wy_cond"], []),
    "RF-008": (["test007_KITEx_gr_ss", "test015_KITEx_wy_ss"], []),
    "RF-009": (["test005_KITEx_gr", "test017_KITEt_gr_opt"], []),
    "RF-011": ([], ["custom_two_sign_convention_regression.py",
                    "custom_two_same_vertex_sign_regression.py",
                    "custom_one_local_reader_regression.py",
                    "custom_output_units_regression.py"]),
    "RF-012": ([], ["custom_two_mod4_regression.py",
                    "custom_output_units_regression.py"]),
}

#: RF-007 and RF-008 each claim a 3D case; these are the Weyl directories.
THREE_DIMENSIONAL = {"RF-007": "test016_KITEx_wy_cond",
                     "RF-008": "test015_KITEx_wy_ss"}

LEGACY = sorted({d for dirs, _ in COVERAGE.values() for d in dirs})


@pytest.mark.parametrize("capability", sorted(COVERAGE))
def test_every_capability_names_something(capability):
    named_tests_exist(capability, *COVERAGE[capability])


@pytest.mark.parametrize("capability", sorted(THREE_DIMENSIONAL))
def test_the_conductivities_are_covered_in_three_dimensions(capability):
    """Read from the stored configuration, not from `wy` meaning Weyl."""
    h5py = pytest.importorskip("h5py")
    name = THREE_DIMENSIONAL[capability]
    with h5py.File(TESTS / name / "configORIG.h5", "r") as f:
        dim = int(f["DIM"][()])
    assert dim == 3, (
        f"{name} is {dim}-dimensional, and it is the only directory carrying "
        f"{capability} in three dimensions")


@pytest.mark.skipif(not BINARIES_BUILT, reason="needs a built KITEx and KITE-tools")
@pytest.mark.parametrize("directory", LEGACY)
def test_the_named_directory_passes(directory):
    assert_the_directory_passes(directory)
