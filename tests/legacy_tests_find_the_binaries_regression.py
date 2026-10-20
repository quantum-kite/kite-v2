"""Every legacy test directory can find the binaries the build produces.

The `tests/test0*` and `tests/large*` directories predate the CMake layout and
called `../KITEx`, i.e. `tests/KITEx`, which nothing has emitted for a long
time. Every one of them failed at the first command, so the whole suite was
dead weight rather than a check.

They now reach the binaries through `$KITEX` and `$KITE_TOOLS`, falling back to
`../../build/` when a test is run on its own. `tests/start_tests.sh` resolves
both to absolute paths before the loop, so a test invoked from elsewhere still
finds them and a missing build reports itself once rather than once per test.

What is asserted here is that every directory resolves to the real binaries --
not that the tests pass. Several do not, for reasons that belong to their own
tasks: stale `configORIG.h5` fixtures that predate the `/Seed0` and `/Seed1`
fields, and directories that assert nothing at all.
"""

import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "build"

#: Directories carrying the shell-driven legacy tests.
LEGACY = sorted([p for p in (REPO / "tests").glob("test0*") if p.is_dir()]
                + [p for p in (REPO / "tests").glob("large*") if p.is_dir()])

#: How a test names a binary: through the variable, with a fallback for when the
#: test is run on its own. The fallback is captured as written, whatever it is,
#: so one pointing somewhere other than the build fails rather than going
#: unrecognised. A bare `../KITEx` is the drift this guards against.
REFERENCE = re.compile(r"\$\{(KITEX|KITE_TOOLS):-([^}]+)\}")
BARE = re.compile(r"(?<!build/)\.\./KITE(?:x|-tools)\b")


def _sources(directory):
    return [p for p in (directory / "test.sh", directory / "test.py") if p.exists()]


def test_there_are_legacy_directories_to_check():
    assert len(LEGACY) > 15, f"only {len(LEGACY)} legacy directories found; the glob is wrong"


@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_no_test_still_points_at_the_old_location(directory):
    for source in _sources(directory):
        text = source.read_text()
        stale = BARE.findall(text)
        assert not stale, (
            f"{source.relative_to(REPO)} still calls {stale[0]}, which would be "
            "tests/KITEx -- the build has never put a binary there")


@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_every_reference_resolves_to_the_build_output(directory):
    found = False
    for source in _sources(directory):
        for _, default in REFERENCE.findall(source.read_text()):
            found = True
            resolved = (directory / default).resolve()
            assert resolved.parent == BUILD.resolve(), (
                f"{source.relative_to(REPO)} resolves to {resolved}, not the build directory")
    if not found:
        pytest.skip(f"{directory.name} invokes neither binary")


@pytest.mark.skipif(not (BUILD / "KITEx").exists(),
                    reason="needs a build to resolve against")
@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_the_resolved_binary_is_executable(directory):
    for source in _sources(directory):
        for _, default in REFERENCE.findall(source.read_text()):
            resolved = (directory / default).resolve()
            assert resolved.exists() and os.access(resolved, os.X_OK), (
                f"{source.relative_to(REPO)} points at {resolved}, which is not executable")


def test_the_driver_refuses_to_run_without_the_binaries():
    """One clear message, instead of every test failing on its own."""
    driver = (REPO / "tests" / "start_tests.sh").read_text()
    assert "KITEX" in driver and "KITE_TOOLS" in driver, (
        "start_tests.sh no longer resolves the binaries for the tests it runs")
    assert "export KITEX" in driver, "the resolved paths are not exported to the tests"
    assert re.search(r"exit 1", driver), (
        "start_tests.sh does not stop when the binaries are missing")
