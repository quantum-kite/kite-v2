"""A flag KITE-tools advertises in `--help` must do something.

`--info` was listed in `--help` as printing information about the compilation
process, was accepted by the parser, and printed nothing: `printInfo()` had its
whole body commented out while the generated header it needed was built and on
the include path. A flag that is advertised and inert teaches the reader to
distrust the rest of the list, which is the real cost.

This checks the class rather than that one flag. Every `--name` the binary's own
`--help` advertises must either produce output when invoked, or be a calculation
keyword the parser recognises in `valid_keys`. A new flag added to the help text
and left unwired fails here.
"""

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TOOLS = REPO / "build" / "KITE-tools"
PARSE = REPO / "tools" / "Src" / "Tools" / "parse_input.cpp"

pytestmark = pytest.mark.skipif(not TOOLS.exists(), reason="needs a built KITE-tools")


def _help():
    return subprocess.run([str(TOOLS), "--help"], capture_output=True, text=True).stdout


def _valid_keys():
    m = re.search(r"valid_keys\s*(?:=|\{)(.*?)\};", PARSE.read_text(), re.S)
    return set(re.findall(r'"(--[A-Za-z][A-Za-z0-9_]*)"', m.group(1))) if m else set()


def _advertised():
    """Flags named in the help text's own option lines, not in its prose."""
    return {m.group(1) for line in _help().splitlines()
            if (m := re.match(r"\s*(--[A-Za-z][A-Za-z0-9_]*)\s", line))}


def test_the_help_text_advertises_something():
    assert _advertised(), "no flags parsed out of --help; the format changed"


@pytest.mark.parametrize("flag", sorted(_advertised()) if TOOLS.exists() else [])
def test_every_advertised_flag_is_wired(flag):
    """Either it is a calculation keyword, or invoking it produces output."""
    if flag in _valid_keys():
        return                       # a calculation selector; needs an HDF5 to do anything
    done = subprocess.run([str(TOOLS), flag], capture_output=True, text=True)
    assert (done.stdout + done.stderr).strip(), (
        f"{flag} is advertised in --help but produced no output at all. It is not "
        "a calculation keyword either, so nothing it claims to do is reachable")


def test_info_reports_the_build():
    """The case that prompted this file: --info must name the build."""
    for form in ("--info", "-i"):
        out = subprocess.run([str(TOOLS), form], capture_output=True, text=True).stdout
        assert "compilation" in out.lower(), f"{form} printed nothing about the build"
        assert out.count("\n") >= 3, (
            f"{form} printed {out.count(chr(10))} lines; expected the machine, the "
            "operating system and the date")
