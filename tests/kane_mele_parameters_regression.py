"""The Kane-Mele numbers in the documentation must match the generator.

docs/documentation/examples/custom_vertex_operators.md quotes t=1, t2=0.15,
64x64, 256 moments, and a plateau region |E| < 3*sqrt(3)*t2 = 0.78|t|. Those
numbers are arithmetic on the generator's defaults, so if the defaults move the
documentation goes stale with nothing to catch it.

Earlier this file also asserted that particular sentences appeared in the
post-processing script's docstring. That pinned wording rather than physics and
is what kept the docstring long, so those checks are gone; the plateau
discussion belongs in the documentation, not in the script.
"""

import ast
import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
GENERATOR = REPO / "examples" / "kane_mele_spin_hall.py"
DOC = REPO / "docs" / "documentation" / "examples" / "custom_vertex_operators.md"


def _generator_defaults():
    tree = ast.parse(GENERATOR.read_text())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    names = [a.arg for a in main.args.args][-len(main.args.defaults):]
    return dict(zip(names, [getattr(d, "value", None) for d in main.args.defaults]))


def test_generator_still_ships_the_documented_parameters():
    defaults = _generator_defaults()
    assert defaults["t2"] == 0.15
    assert defaults["moments"] == 256


def test_documented_plateau_width_follows_from_that_t2():
    """The doc's 0.78|t| is 3*sqrt(3)*t2, not an independent measurement."""
    t2 = _generator_defaults()["t2"]
    assert 3 * np.sqrt(3) * t2 == pytest.approx(0.78, abs=5e-3)
    assert "0.78" in DOC.read_text(), "the documented plateau width moved"


def test_documentation_quotes_the_parameters_the_generator_uses():
    text = DOC.read_text()
    defaults = _generator_defaults()
    assert re.search(r"t_2\s*=\s*0\.15", text), "t2 in the docs no longer matches"
    assert str(defaults["moments"]) in text, "the moment count in the docs moved"
