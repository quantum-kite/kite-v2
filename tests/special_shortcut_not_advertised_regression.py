"""The second-order shortcut stays in the code, and out of the documentation.

`conductivity_optical_nonlinear` accepts a `special` argument that selects an
internal path computing two of the four Chebyshev objects instead of all four.
It is a materials-specific simplification, not a numerical one -- the C++
comment calls it "nonlinear but only has simple objects that need calculating"
-- so it is not something a user can choose on general grounds, and a reader who
meets it in an API table has no way to judge when it applies.

It keeps working, because scripts in the wild pass it. What it must not do is
appear anywhere that invites a user to set it (AC-NLO-02).

The check is a search rather than a fixture, because the failure mode is
someone reintroducing the mention somewhere new -- a new example, a new page --
not an edit to a known line. The package's own implementation is exempt by
construction: for Python it inspects docstrings only, never the code that
honours the argument.
"""

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: The argument shown to a reader: as a keyword in a call (`special=1`), or as
#: a documented parameter in the API tables, which mark identifiers `#!python x`.
#: Matching only `special=` would miss the table row, which carries no `=`.
_ADVERTISED = re.compile(r"special\s*=|#!python\s+special\b|`special`")


def _pages():
    """Markdown anywhere a user reads, and every example script in full."""
    for root in (REPO / "docs", REPO / "examples"):
        for path in sorted(root.rglob("*")):
            if path.suffix in {".md", ".rst", ".py"}:
                yield path


def _docstrings(path):
    """Every docstring in a Python module, with its line number."""
    tree = ast.parse(path.read_text(encoding="utf-8", errors="ignore"))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            text = ast.get_docstring(node)
            if text:
                yield getattr(node, "lineno", 1), text


def test_there_are_pages_to_search():
    """A glob that matched nothing would make the next tests vacuous."""
    assert sum(1 for _ in _pages()) > 50


@pytest.mark.parametrize("path", [p.relative_to(REPO).as_posix() for p in _pages()])
def test_no_page_or_example_advertises_the_shortcut(path):
    for line_no, line in enumerate((REPO / path).read_text(
            encoding="utf-8", errors="ignore").splitlines(), 1):
        if _ADVERTISED.search(line):
            pytest.fail(
                f"{path}:{line_no} shows `special=` to the user:\n    {line.strip()}\n"
                "The argument still works, but it selects a materials-specific "
                "path whose applicability a reader cannot judge, so it is not "
                "documented (AC-NLO-02)")


def test_the_api_docstring_does_not_mention_the_shortcut():
    """The package may use `special`; it may not explain it to the caller."""
    module = REPO / "python" / "kite" / "__init__.py"
    offenders = [(n, t) for n, t in _docstrings(module) if "special" in t.lower()]
    assert not offenders, (
        "a docstring in python/kite/__init__.py mentions the shortcut: "
        + "; ".join(f"line {n}" for n, _ in offenders))
