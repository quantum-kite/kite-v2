"""Every cross-reference in the documentation must resolve.

Three failure modes, none of which `mkdocs --strict` catches on its own:

* a reference-style link `[text][ref]` whose `[ref]:` definition is missing.
  Markdown renders that as literal text -- the reader sees `[text][ref]` on the
  page and the build stays green.
* a relative link to a `.md` file that no longer exists.
* a link to a personal fork rather than to this repository, which sends readers
  to a tree nobody maintains. `docs/documentation/examples/spectral_function.md`
  pointed at `tgrappoport/kite` for a script that ships here.

Links out to genuinely external projects are fine and are listed below.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DOCS = REPO / "docs"

#: Upstream projects the docs legitimately point at.
EXTERNAL_OWNERS = {"quantum-kite", "dean0x7d", "conda-forge"}

PAGES = sorted(DOCS.rglob("*.md"))


def _body(text):
    """Strip fenced code; a `[a][b]` inside an example is not a link."""
    return re.sub(r"```.*?```", "", text, flags=re.S)


@pytest.mark.parametrize("page", PAGES, ids=lambda p: str(p.relative_to(DOCS)))
def test_every_reference_style_link_is_defined(page):
    text = page.read_text(errors="ignore")
    defined = set(re.findall(r"^\[([^\]]+)\]:\s*\S", text, re.M))
    used = {r for r in re.findall(r"\]\[([^\]]+)\]", _body(text)) if r}
    missing = sorted(used - defined)
    assert not missing, (
        f"{page.relative_to(DOCS)} uses {missing} with no matching definition; "
        "the page renders the raw bracket text and the build still passes")


@pytest.mark.parametrize("page", PAGES, ids=lambda p: str(p.relative_to(DOCS)))
def test_every_relative_markdown_target_exists(page):
    broken = []
    for target in re.findall(r"\]\(([^)#\s]+\.md)(?:#[^)]*)?\)", _body(page.read_text(errors="ignore"))):
        if target.startswith(("http", "mailto")):
            continue
        if not (page.parent / target).resolve().exists():
            broken.append(target)
    assert not broken, f"{page.relative_to(DOCS)} links to missing pages: {broken}"


def test_no_page_links_to_a_personal_fork():
    offenders = []
    for page in PAGES:
        for owner in re.findall(r"https?://github\.com/([\w\-]+)/[\w\-]+", page.read_text(errors="ignore")):
            if owner not in EXTERNAL_OWNERS:
                offenders.append((str(page.relative_to(DOCS)), owner))
    assert not offenders, (
        f"the documentation links to non-upstream GitHub owners {offenders}. A "
        "link to a personal fork sends readers to a tree this project does not "
        "maintain")
