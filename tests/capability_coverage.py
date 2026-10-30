"""Shared machinery for the capability-coverage files.

Three of them say which test discharges which requirement and then run it, so
the part that runs a legacy directory lives here rather than three times over.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTS = REPO / "tests"
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

BINARIES_BUILT = KITEX.exists() and KITE_TOOLS.exists()


def named_tests_exist(capability, legacy, modern):
    """Raise AssertionError naming the capability if anything it points at is gone."""
    assert legacy or modern, f"{capability} names no test at all"
    for name in legacy:
        assert (TESTS / name / "test.py").exists(), (
            f"{capability} names {name}, which has no test.py")
    for name in modern:
        assert (TESTS / name).exists(), (
            f"{capability} names {name}, which is not in tests/")


def run_legacy_directory(directory):
    """The quick path: restore the stored input, run, return stdout.

    `config.h5` is gitignored while `configORIG.h5` and `configREF.h5` are
    tracked, so running in place leaves the working tree clean.
    """
    d = TESTS / directory
    shutil.copy(d / "configORIG.h5", d / "config.h5")
    env = dict(os.environ, KITEX=str(KITEX), KITE_TOOLS=str(KITE_TOOLS))
    proc = subprocess.run([sys.executable, "test.py"], cwd=d, env=env,
                          capture_output=True, text=True)
    assert proc.returncode == 0, f"{directory}/test.py exited {proc.returncode}"
    return proc.stdout


def assert_the_directory_passes(directory):
    out = run_legacy_directory(directory)
    assert "ERROR" not in out, (
        f"{directory} could not run KITEx: {out.strip()!r}")
    assert "Problem" not in out, (
        f"{directory} disagrees with its stored reference: {out.strip()!r}")
    assert "OK" in out, (
        f"{directory} printed no verdict at all: {out.strip()!r}. A test that "
        "says nothing is indistinguishable from one that passed")
