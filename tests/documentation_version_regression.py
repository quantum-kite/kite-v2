"""Every place that declares a KITE version must declare the same one.

The documentation states which release it describes, in four places that have
no mechanical link to each other: `VERSION.md`, `README.md`,
`docs/background/index.md` and `mkdocs.yml`. They agree today. They will not
stay agreeing through a release unless something checks, and a page left
claiming the previous version is exactly the failure this guards -- a reader
cannot tell which KITE the page is about, and the capability list that follows
it is the part they will act on.

This deliberately does not assert *which* version. Bumping is a release
decision, not a defect, so the test only requires that the four agree.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

#: (file, regex capturing the version it declares)
DECLARATIONS = [
    ("VERSION.md", r"Version\s+(\d+\.\d+)"),
    ("README.md", r"KITE v(\d+\.\d+) \(this version\)"),
    ("docs/background/index.md", r"latest release \(version (\d+\.\d+)\)"),
    ("mkdocs.yml", r'^\s*version:\s*"(\d+\.\d+)"'),
]


def _declared(path, pattern):
    text = (REPO / path).read_text()
    match = re.search(pattern, text, re.M)
    assert match, (
        f"{path} no longer declares a version matching {pattern!r}; the "
        "documentation would then not say which release it describes")
    return match.group(1)


def test_every_declaration_is_present_and_parseable():
    for path, pattern in DECLARATIONS:
        _declared(path, pattern)


def test_all_the_declarations_agree():
    found = {path: _declared(path, pattern) for path, pattern in DECLARATIONS}
    distinct = set(found.values())
    assert len(distinct) == 1, (
        "the documented version differs between files: "
        + ", ".join(f"{p} says {v}" for p, v in sorted(found.items()))
        + ". A page left on the previous version misdescribes the capabilities "
          "listed beneath it")


def test_the_site_config_carries_the_version():
    """So the marker travels with the built site, not only with the prose."""
    text = (REPO / "mkdocs.yml").read_text()
    assert re.search(r'^\s*version:\s*"\d+\.\d+"', text, re.M), (
        "mkdocs.yml no longer records the documented release under `extra:`")
