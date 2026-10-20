"""Every shipped example must appear in the catalogue.

`examples/README.md` is the catalogue, and it is also a documentation page --
`docs/documentation/more_examples/additional_examples.md` is a symlink to it, so
one file serves both the repository and the site. A script added to the folder
without an entry is invisible to a reader browsing the documentation, and a
script removed without its entry going too leaves the page describing something
that is not there.

Not every entry is a full write-up. A handful of scripts are shared helpers or
test fixtures rather than examples, and the catalogue says so in prose instead
of giving them a table row; that is a documented decision, not an omission. What
this checks is that nothing is simply absent.
"""

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"
CATALOGUE = EXAMPLES / "README.md"

#: Library modules, not scripts: imported by the examples, never run alone.
HELPERS = {"myaux.py", "kite_style.py"}


def _scripts():
    return sorted(p.name for p in EXAMPLES.glob("*.py") if p.name not in HELPERS)


def test_the_catalogue_is_the_documentation_page():
    """One file, two readers. If the symlink breaks they drift apart silently."""
    page = REPO / "docs" / "documentation" / "more_examples" / "additional_examples.md"
    assert page.is_symlink(), (
        f"{page.relative_to(REPO)} is no longer a symlink to the catalogue; the "
        "repository and the site can now disagree about what ships")
    assert page.resolve() == CATALOGUE.resolve(), (
        f"the documentation page points at {page.resolve()}, not the catalogue")


@pytest.mark.parametrize("script", _scripts())
def test_every_example_appears_in_the_catalogue(script):
    text = CATALOGUE.read_text()
    assert script in text, (
        f"examples/{script} ships but is not mentioned in examples/README.md, so "
        "it does not appear in the documentation either")


def test_the_catalogue_mentions_no_script_that_was_removed():
    """The other direction: a deleted example must not linger in the page.

    Every `.py` under `examples/` counts, not just the top level, because the
    page also describes the `pybinding` and `paper` trees. Matching on the bare
    filename is deliberate -- the page names scripts, not paths.
    """
    import re

    shipped = {p.name for p in EXAMPLES.rglob("*.py")}
    named = set(re.findall(r"`([A-Za-z0-9_]+\.py)`", CATALOGUE.read_text()))

    stale = sorted(named - shipped)
    assert not stale, (
        f"the catalogue describes scripts that are no longer shipped: {stale}")
