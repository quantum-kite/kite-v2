"""The lattice repository, built and checked material by material.

graphene, phosphorene and the group-6 transition-metal dichalcogenides ship as
a repository that the examples build on, and nothing tested it directly: a
change that broke a material would surface eventually as a broken example,
which is a slow and indirect way to find out. This file builds every lattice
the repository offers and checks its shape against what the parameter tables
say it should be.

The TMD check walks all six materials rather than one, because the parameters
are a per-material table and a typo in one row is exactly the kind of thing
that builds cleanly and is wrong.
"""

import pytest

#: material -> the 3-band model puts three orbitals in the cell for each
TMD_MATERIALS = ["MoS2", "WS2", "MoSe2", "WSe2", "MoTe2", "WTe2"]


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
