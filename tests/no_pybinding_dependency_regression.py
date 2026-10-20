"""Nothing that ships as a test or an example requires pybinding.

`pybinding` is an optional extra (`pyproject.toml`'s `[project.optional-dependencies]`),
and it is not installed in the environment `environment.yml` describes. Anything
that imports it is therefore dead for most users: the legacy test directories
could not regenerate their inputs, and one of them shipped no input at all and so
could not run.

Everything those files used has a counterpart KITE ships itself -- `pb.Lattice`
is `kite.lattice.Lattice` with the same constructor and the same
`add_sublattices`/`add_hoppings`, `pb.results.make_path` is
`kite.visualize.make_path` (which returns a tuple, so callers index `[0]`), and
`pybinding.repository.graphene`'s `a`/`a_cc` are in `kite.repository.graphene`
with identical values.

`examples/pybinding/` is deliberately exempt: it exists to show KITE driven from
pybinding, so importing it is the point. `examples/paper/` is reproduction
material for a published figure set and is covered by its own task.
"""

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: The one tree where importing pybinding is the purpose, not a dependency.
EXEMPT = (REPO / "examples" / "pybinding", REPO / "examples" / "paper")


def _python_files():
    for root in (REPO / "tests", REPO / "tools" / "tests", REPO / "examples"):
        for path in root.rglob("*.py"):
            if not any(str(path).startswith(str(e)) for e in EXEMPT):
                yield path


def _imports_pybinding(path):
    """A real import, not a mention in a docstring or a comment."""
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(a.name.split(".")[0] == "pybinding" for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] == "pybinding":
                return True
    return False


FILES = sorted(_python_files())


def test_there_are_files_to_check():
    assert len(FILES) > 60, f"only {len(FILES)} files found; the walk is wrong"


@pytest.mark.parametrize("path", FILES, ids=lambda p: str(p.relative_to(REPO)))
def test_no_shipped_test_or_example_imports_pybinding(path):
    assert not _imports_pybinding(path), (
        f"{path.relative_to(REPO)} imports pybinding, which is an optional extra "
        "and is absent from the documented environment")


def test_the_replacements_exist_and_agree_with_what_they_replaced():
    """If these drift, the ported configs build a different lattice in silence."""
    from kite import lattice as latt, visualize
    from kite.repository import graphene
    import inspect

    params = inspect.signature(latt.Lattice.__init__).parameters
    for name in ("a1", "a2", "a3"):
        assert name in params, f"kite.lattice.Lattice lost the {name} argument"
    for method in ("add_sublattices", "add_hoppings"):
        assert hasattr(latt.Lattice, method), f"kite.lattice.Lattice lost {method}"

    assert callable(visualize.make_path)
    # pybinding's graphene constants, which test013 used directly.
    assert graphene.a == pytest.approx(0.24595)
    assert graphene.a_cc == pytest.approx(0.142)


@pytest.mark.parametrize("name", ["test001_KITEx_sq", "test004_KITEx_cb", "test013_KITEx_km_arpes"])
def test_a_ported_config_still_builds_its_lattice(name, tmp_path):
    """One 2D, one 3D, one that used pybinding's graphene constants."""
    import shutil
    import subprocess
    import sys
    import os

    source = REPO / "tests" / name
    for script in source.glob("*.py"):
        shutil.copy(script, tmp_path / script.name)
    env = {**os.environ, "PYTHONPATH": str(REPO / "python"), "MPLBACKEND": "Agg"}
    done = subprocess.run([sys.executable, "config.py"], cwd=tmp_path, env=env,
                          capture_output=True, text=True)
    assert done.returncode == 0, f"{name}/config.py failed:\n{done.stderr[-600:]}"
    assert (tmp_path / "config.h5").exists(), f"{name}/config.py wrote no config.h5"
