"""A structural-disorder value is an increment, and the page must say so.

`add_structural_disorder` takes a number per bond and per site. The number is
**added** to the clean lattice, not substituted for it: `KPM_VectorBasis.cpp:136`
accumulates `phi0[k1] += new_hopping * phiM1[k2]` on top of the contribution the
regular hopping loop has already made.

That is the design -- you pass the change -- but it is easy to read the other
way, and a shipped example once wrote `timp = -0.0` expecting to sever a bond.
It severed nothing. To cancel a hopping of `-t` you pass `+t`.

A behavioural test is not viable here: structural defects are placed at random
by concentration, and the resulting DOS difference is smaller than the
stochastic spread. I measured both -- clean-vs-clean across two seeds gives
6.2e-02, while clean-vs-(value 0.0) gives 4.6e-02 -- so a run-based check would
be reading noise. This pins the source instead.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_the_defect_hopping_is_accumulated_not_assigned():
    """`+=`, not `=`. The whole semantics rests on this operator."""
    source = (REPO / "Src" / "Vector" / "KPM_VectorBasis.cpp").read_text()
    assert re.search(r"phi0\[k1\]\s*\+=\s*vel \* tmp", source), (
        "the structural-defect contribution is no longer accumulated with += "
        "into phi0. If it is assigned instead, the value became a replacement "
        "and docs/documentation/disorder.md is wrong")

    header = (REPO / "Src" / "Hamiltonian" / "HamiltonianDefects.hpp").read_text()
    assert header.count("phi0[k1] +=") >= 2, (
        "HamiltonianDefects.hpp no longer accumulates the defect hopping; the "
        "velocity and non-velocity branches both used to")


def test_the_page_states_the_increment_semantics_and_the_units():
    page = (REPO / "docs" / "documentation" / "disorder.md").read_text()
    assert "increments added to the clean lattice" in page, (
        "disorder.md no longer says the structural-disorder values are "
        "increments; a reader will pass the intended new hopping instead of "
        "the change, as a shipped example once did")
    assert "same units as the hoppings" in page, (
        "disorder.md no longer states the units of the structural-disorder "
        "values")
    assert "a value of $0$ changes nothing" in page, (
        "disorder.md no longer gives the consequence that makes the semantics "
        "concrete")
