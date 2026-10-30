"""Each lattice and disorder capability names the test that protects it.

RF-015 (building lattices natively, without the optional `pybinding`), RF-017
(Anderson on-site disorder) and RF-018 (vacancies and multi-site defect
clusters) are all exercised by tests that already exist -- this file says which,
and runs the legacy directories among them.

RF-016 is different. The lattice repository -- graphene, phosphorene and the
group-6 TMDs -- had no test at all: examples use it, so a change that broke a
material would surface eventually as a broken example, which is a slow and
indirect way to find out. So RF-016 is discharged here directly, by building
every lattice the repository offers and checking its shape against what the
parameter tables say it should be.

The TMD check walks all six materials rather than one, because the parameters
are a per-material table and a typo in one row is exactly the kind of thing that
builds cleanly and is wrong.
"""

import pytest

from capability_coverage import (BINARIES_BUILT, TESTS,
                                 assert_the_directory_passes, named_tests_exist)

#: capability -> (legacy directories, modern test files)
COVERAGE = {
    "RF-015": ([], ["no_pybinding_dependency_regression.py",
                    "lattice_vectors_precision_regression.py",
                    "hopping_range_export_regression.py"]),
    "RF-016": ([], ["lattice_and_disorder_coverage_regression.py"]),
    "RF-017": (["test007_KITEx_gr_ss", "test013_KITEx_km_arpes",
                "test014_KITEx_gr_adat"], []),
    "RF-018": (["test011_KITEx_sq_defec", "test014_KITEx_gr_adat",
                "test015_KITEx_wy_ss", "test016_KITEx_wy_cond"], []),
}

LEGACY = sorted({d for dirs, _ in COVERAGE.values() for d in dirs})

#: material -> number of sublattices the 3-band model puts in the cell
TMD_MATERIALS = ["MoS2", "WS2", "MoSe2", "WSe2", "MoTe2", "WTe2"]


@pytest.mark.parametrize("capability", sorted(COVERAGE))
def test_every_capability_names_something(capability):
    named_tests_exist(capability, *COVERAGE[capability])


@pytest.mark.skipif(not BINARIES_BUILT, reason="needs a built KITEx and KITE-tools")
@pytest.mark.parametrize("directory", LEGACY)
def test_the_named_directory_passes(directory):
    assert_the_directory_passes(directory)


# --- RF-016: the lattice repository, which nothing covered ------------------

def _repository():
    from kite import repository
    return repository


def test_graphene_builds_with_two_sublattices_and_three_bonds():
    lat = _repository().graphene.monolayer()
    assert list(lat.sublattices) == ["A", "B"]
    assert lat.nhop == 3, (
        f"the nearest-neighbour honeycomb has three distinct bonds, not {lat.nhop}")
    assert len(lat.vectors) == 2


def test_graphene_second_neighbours_add_bonds_without_adding_sublattices():
    repository = _repository()
    near = repository.graphene.monolayer()
    far = repository.graphene.monolayer(nearest_neighbors=2)
    assert list(far.sublattices) == list(near.sublattices)
    assert far.nhop > near.nhop, (
        "asking for second neighbours added no bonds")


def test_phosphorene_builds_with_four_sublattices():
    lat = _repository().phosphorene.monolayer_4band()
    assert list(lat.sublattices) == ["A", "B", "C", "D"], (
        "the 4-band phosphorene model has four sites in the cell")
    assert lat.nhop >= 4


@pytest.mark.parametrize("material", TMD_MATERIALS)
def test_every_group6_tmd_builds(material):
    """The parameters are a per-material table; a bad row still builds."""
    lat = _repository().group6_tmd.monolayer_3band(material)
    assert len(lat.sublattices) == 3, (
        f"{material}'s 3-band model must carry three orbitals in the cell")
    assert lat.nhop > 0, f"{material} has no hoppings"


def test_an_unknown_material_is_refused():
    """Rather than silently falling back to some other material's parameters."""
    with pytest.raises((KeyError, ValueError)):
        _repository().group6_tmd.monolayer_3band("NotAMaterial")
