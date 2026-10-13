"""No shipped model may declare a hopping beyond the first neighbouring cell.

KITE encodes a hopping's relative cell index in a fixed-width address: each
component is stored in one base-3 digit, so the only representable values are
-1, 0 and +1. A model declaring (-2, +1) does not fail -- it is silently
encoded as something else, and KITEx computes a different Hamiltonian than the
script describes.

`tests/hopping_range_export_regression.py` guards the export: it checks that
`config_system` refuses such a model. This file is the other half, and the one
the task asks for: it sweeps the models actually shipped, so a lattice that
declares an out-of-range index is caught even if nobody exports it in CI.

Every example is imported and its lattice captured, rather than grepping for
literals, because a relative index can be built in a loop or read from a variable.

`examples/paper/` is excluded: it is reproduction material for a published
figure set, kept as it was rather than maintained alongside the examples, and
it needs `pybinding`, which is not a project dependency.
"""

import contextlib
import importlib.util
import io
import sys
from pathlib import Path

import re

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"

#: Base-3 addressing: one digit per dimension, so only -1, 0, +1 are representable.
LIMIT = 1

SKIP = {"myaux.py", "kite_style.py"}


class _Captured(Exception):
    """Raised to stop an example once its lattice is in hand."""


def _lattices():
    """Every example's lattice, captured at the config_system call.

    Nothing is written to disk and no calculation runs: the spy raises as soon
    as the lattice is available.
    """
    import kite

    # Importing arbitrary example modules leaves their imports behind in
    # sys.modules, which is visible to every test that runs afterwards --
    # tests/estimate_bounds_regression.py asserts that `pybinding` was never
    # imported, and this sweep reaching a paper script that imports it made
    # that test fail. Snapshot and restore rather than leaving the pollution.
    before = set(sys.modules)

    found = []
    for path in sorted(EXAMPLES.rglob("*.py")):
        if path.name in SKIP or {"pybinding", "paper"} & set(path.parts):
            continue
        captured = {}
        real = kite.config_system

        def spy(lattice, *args, **kwargs):
            captured["lattice"] = lattice
            raise _Captured()

        kite.config_system = spy
        sys.path.insert(0, str(path.parent))
        try:
            spec = importlib.util.spec_from_file_location("swept_" + path.stem, path)
            module = importlib.util.module_from_spec(spec)
            with contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()):
                try:
                    spec.loader.exec_module(module)
                except _Captured:
                    pass
                except BaseException:
                    continue
                if "lattice" not in captured and hasattr(module, "main"):
                    try:
                        module.main()
                    except _Captured:
                        pass
                    except BaseException:
                        continue
        finally:
            kite.config_system = real
            sys.path.remove(str(path.parent))
        if "lattice" in captured:
            found.append((str(path.relative_to(REPO)), captured["lattice"]))

    for name in set(sys.modules) - before:
        del sys.modules[name]
    return found


def _out_of_range(lattice):
    """Relative indices with any component outside [-1, +1]."""
    offenders = []
    for family in getattr(lattice, "hoppings", {}).values():
        for term in getattr(family, "terms", []):
            index = np.asarray(term.relative_index, dtype=int).ravel()
            if np.any(np.abs(index) > LIMIT):
                offenders.append(tuple(int(i) for i in index))
    return sorted(set(offenders))


MODELS = _lattices()


def test_the_sweep_captured_the_shipped_models():
    """A sweep that silently captures nothing would pass forever."""
    assert len(MODELS) >= 20, (
        f"only {len(MODELS)} lattices captured; the sweep is not working")


@pytest.mark.parametrize("name,lattice", MODELS, ids=[n for n, _ in MODELS])
def test_no_shipped_model_reaches_beyond_the_first_neighbour_cell(name, lattice):
    offenders = _out_of_range(lattice)
    assert not offenders, (
        f"{name} declares relative indices {offenders}, outside [-1, +1]. KITE "
        "stores each component in one base-3 digit, so this is encoded as a "
        "different hopping and the Hamiltonian computed is not the one declared.")


def test_the_sweep_would_catch_an_out_of_range_model():
    """The check itself, against a lattice built to violate it."""
    from kite import lattice as latt

    a_cc = 1 / np.sqrt(3)
    lat = latt.Lattice(a1=[1, 0], a2=[0.5, np.sqrt(3) / 2])
    lat.add_sublattices(("A", [0, -a_cc / 2], 0.0), ("B", [0, a_cc / 2], 0.0))
    lat.add_hoppings(([0, 0], "A", "B", -1.0), ([-2, 1], "A", "B", -0.1))
    assert _out_of_range(lat) == [(-2, 1)]
