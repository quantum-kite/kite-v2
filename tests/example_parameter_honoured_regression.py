"""Regression: an example must not silently override its own parameters.

`dos_vacancies.honeycomb_lattice(onsite, t)` declared `t=1`, then reassigned
`t = 2.8` unconditionally two lines later. Three statements disagreed: the
signature said 1, the body said 2.8, and the module docstring said "Units:
Energy in units of hopping, |t| = 1".

The consequence was not cosmetic. The example uses automatic scaling, so the
exported `EnergyScale` is derived from the hoppings -- it was 16.8 while every
energy a reader computed from the docstring was in units of t = 1. Any acceptance
criterion written against the documented units would have been wrong by 2.8.

This test asserts the caller's value reaches the exported Hamiltonian, which is
the property that was broken, rather than asserting the literal 1 -- so it keeps
working if the example's default ever changes deliberately.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "dos_vacancies.py"


def _load_example():
    sys.path.insert(0, str(REPO / "python"))
    spec = importlib.util.spec_from_file_location("dos_vacancies", EXAMPLE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("t", [1.0, 2.8, 0.5])
def test_hopping_parameter_reaches_the_lattice(t):
    """Whatever the caller passes is the hopping the lattice is built with."""
    mod = _load_example()
    lattice = mod.honeycomb_lattice(t=t)
    # HoppingFamily.energy is a (1, 1) complex array, not a scalar, and the
    # honeycomb hopping is stored as -t.
    magnitudes = sorted({
        float(np.abs(np.asarray(hf.energy)).max())
        for hf in lattice.hoppings.values()
    })
    assert magnitudes == pytest.approx([t]), (
        f"caller passed t={t} but the lattice was built with {magnitudes}"
    )


def test_docstring_unit_claim_matches_the_default():
    """The module documents |t| = 1; the default must actually be 1."""
    import inspect

    mod = _load_example()
    sig = inspect.signature(mod.honeycomb_lattice)
    assert sig.parameters["t"].default == 1
    assert "|t| = 1" in mod.__doc__

    src = inspect.getsource(mod.honeycomb_lattice)
    assert "t = 2.8" not in src, "the parameter is being overridden in the body again"
