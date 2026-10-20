"""No legacy test directory reports success without asserting something.

A directory that runs, prints nothing alarming and exits is indistinguishable
from one that passed, which is how `test005`'s zero-against-zero comparison and
`test017`'s commented-out assertion survived for years. `test000_speed` was a
stopwatch whose pass criterion lived in a comment; `test_12` held a description
and no test at all. Both are gone.

What this pins is structural: every directory the driver runs has a `test.py`,
that `test.py` can print both verdicts, and `test.sh`'s quick path actually
calls it. It does not check that the assertion is *meaningful* -- several
comparisons are still against arrays that are identically zero under
`SEED=ones`, which is recorded separately and left to a later pass.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TESTS = REPO / "tests"

#: The directories tests/start_tests.sh loops over.
LEGACY = sorted(p for p in TESTS.glob("test0*") if p.is_dir())


def test_the_driver_still_loops_over_these():
    driver = (TESTS / "start_tests.sh").read_text()
    assert re.search(r"for\s+i\s+in\s+test0\*", driver), (
        "start_tests.sh no longer loops over test0*; this guard is watching the "
        "wrong set")
    assert len(LEGACY) > 15, f"only {len(LEGACY)} directories found"


@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_the_directory_has_a_test(directory):
    assert (directory / "test.py").exists(), (
        f"{directory.name} has no test.py, so running it can only ever report "
        "success")


@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_the_test_can_report_both_verdicts(directory):
    text = (directory / "test.py").read_text()
    # Matched loosely: several of these print "OK " with a trailing space, or
    # "Problem" followed by the comparison numbers.
    for verdict in ("OK", "Problem"):
        assert re.search(rf"print\([^)]*{verdict}", text), (
            f"{directory.name}/test.py never prints {verdict!r}, so it cannot "
            "distinguish a pass from a failure")


@pytest.mark.parametrize("directory", LEGACY, ids=lambda p: p.name)
def test_the_quick_path_runs_the_test(directory):
    """A test.py that the driver never calls is the same as no test."""
    script = (directory / "test.sh").read_text()
    quick = script.split('"quick"', 1)
    assert len(quick) == 2, f"{directory.name}/test.sh has no quick path"
    body = quick[1]
    calls = [line for line in body.splitlines()
             if "test.py" in line and not line.strip().startswith("#")]
    assert calls, (
        f"{directory.name}/test.sh's quick path does not call test.py "
        "(commented out, or absent)")
