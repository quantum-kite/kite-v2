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


#: Every calculation a user can request, as `Calculation` exposes it. The
#: capability lists in README.md and docs/background/index.md describe these in
#: prose, so the prose and this set drift apart silently -- which is how the
#: lists came to omit five of them.
ENTRY_POINTS = {
    "dos", "ldos", "ldos_map", "spectral_map", "arpes",
    "gaussian_wave_packet", "conductivity_dc", "conductivity_optical",
    "conductivity_optical_nonlinear", "singleshot_conductivity_dc",
    "custom_one", "custom_one_local", "custom_two", "custom_two_local",
    "custom_singleshot_two",
}

#: Methods of `Calculation` that are not advertised capabilities:
#: `add_orbital_*` build a custom operator rather than request a calculation,
#: and `localized_wave_packet` is descoped this release
#: (see `descoped_capability_regression.py`).
NOT_ADVERTISED = {"__init__", "add_orbital_coupling", "add_orbital_index",
                  "localized_wave_packet"}


def test_the_capability_lists_still_describe_every_entry_point():
    """Not a wording check -- a prompt to revisit the prose when the set grows.

    There is no mechanical link between a method name and the sentence that
    advertises it, so this pins the set instead. A new entry point fails here,
    and whoever adds it decides how the two lists should say so.
    """
    source = (REPO / "python" / "kite" / "__init__.py").read_text()
    body = source[source.index("class Calculation:"):source.index("class Configuration:")]
    found = set(re.findall(r"^    def (\w+)\(self,", body, re.M))
    found = {name for name in found if not name.startswith("_")}

    added = sorted(found - ENTRY_POINTS - NOT_ADVERTISED)
    removed = sorted(ENTRY_POINTS - found)
    assert not added, (
        f"Calculation gained {added}. README.md and docs/background/index.md "
        "list what KITE can do; add it there, then add it here")
    assert not removed, (
        f"Calculation no longer offers {removed}, but the capability lists in "
        "README.md and docs/background/index.md still advertise it")


def test_the_site_config_carries_the_version():
    """So the marker travels with the built site, not only with the prose."""
    text = (REPO / "mkdocs.yml").read_text()
    assert re.search(r'^\s*version:\s*"\d+\.\d+"', text, re.M), (
        "mkdocs.yml no longer records the documented release under `extra:`")
