"""Configuring without HDF5 must fail, and say what to do; with it, succeed.

Everything KITE reads and writes travels through HDF5 -- the configuration
file, the Chebyshev moments, the results -- and `"H5Cpp.h"` is included
unconditionally, so a build without it cannot compile. `find_package(HDF5 ...)`
without `REQUIRED` used to report a missing dependency as one status line and
configure successfully anyway; the failure then arrived as compiler errors about
a missing header, long after the line explaining it had scrolled past.

The C++ and high-level components matter as much as the library: many
distributions ship HDF5 without them, and that is the common way to have HDF5
installed and still be unable to build.

This mirrors tests/openmp_required_regression.py, including the shape that file
had to be corrected into: both directions are covered, and whether HDF5 is
available at all is established with a throwaway project rather than with
KITE's own CMakeLists, so a guard broken to fire unconditionally cannot
masquerade as a machine without HDF5.

No minimum version is asserted anywhere -- which version KITE actually needs has
not been established, and inventing a number would exclude working setups.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

needs_cmake = pytest.mark.skipif(shutil.which("cmake") is None,
                                 reason="cmake not on PATH")

PROBE = """cmake_minimum_required(VERSION 3.9)
project(probe C CXX)
find_package(HDF5 COMPONENTS CXX C HL)
if(HDF5_FOUND)
  message(STATUS "PROBE_HDF5_YES")
endif()
"""


def _configure(source, build, *extra):
    return subprocess.run(["cmake", "-S", str(source), "-B", str(build), *extra],
                          capture_output=True, text=True)


@pytest.fixture(scope="module")
def without_hdf5(tmp_path_factory):
    """One configure, shared: it costs seconds and both assertions read it."""
    build = tmp_path_factory.mktemp("no-hdf5")
    return _configure(REPO, build, "-DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON")


@pytest.fixture(scope="module")
def hdf5_is_available(tmp_path_factory):
    """Does this toolchain have HDF5 with CXX/C/HL at all?

    Answered with a minimal project of our own. Asking KITE's CMakeLists would
    make a guard that always fires indistinguishable from a machine with no
    HDF5, and the positive test below would skip itself into uselessness.
    """
    root = tmp_path_factory.mktemp("probe")
    (root / "CMakeLists.txt").write_text(PROBE)
    proc = _configure(root, root / "build")
    return "PROBE_HDF5_YES" in proc.stdout


#: The first line of the FATAL_ERROR, used to tell "this guard fired" apart
#: from "the configure failed somewhere else".
GUARD_TEXT = "KITE requires HDF5"


@needs_cmake
def test_configure_fails_without_hdf5_and_says_what_to_do(without_hdf5):
    """It must fail, and the failure must name the cause and a route out."""
    assert without_hdf5.returncode != 0, (
        "configuration succeeded without HDF5; the build would then fail on a "
        "missing H5Cpp.h instead:\n" + without_hdf5.stdout)
    message = without_hdf5.stdout + without_hdf5.stderr
    assert GUARD_TEXT in message, (
        "the configure failed, but not at the HDF5 guard -- something else "
        f"stopped it first, so this proves nothing about HDF5:\n{message}")
    for expected in ("C++", "environment.yml", "docs/installation.md"):
        assert expected in message, (
            f"the failure does not mention {expected!r}, so a user cannot act on it:\n"
            + message)


@needs_cmake
def test_the_guard_does_not_fire_when_hdf5_is_present(hdf5_is_available, tmp_path):
    """The other direction: the guard must not fire when it should not.

    Without this, a guard rewritten to fire unconditionally -- which would stop
    every build on every platform -- leaves the suite green.

    Asserts the guard's own message is absent rather than that the configure
    succeeded outright. HDF5 is checked at CMakeLists.txt:25 and OpenMP at :87,
    so on a machine whose default compiler has no OpenMP the configure fails
    *after* passing this guard. Requiring exit 0 here would make this test
    report on the OpenMP guard instead of its own.
    """
    if not hdf5_is_available:
        pytest.skip("this toolchain has no HDF5 with CXX/C/HL, so the positive "
                    "path cannot be tested")
    proc = _configure(REPO, tmp_path / "with-hdf5")
    message = proc.stdout + proc.stderr
    assert GUARD_TEXT not in message, (
        "HDF5 is available on this toolchain, but the guard fired anyway -- it "
        f"would stop every build:\n{message}")


def test_the_cxx_header_is_included_unconditionally():
    """The guard's premise, read from the source.

    Failing hard is only justified while no translation unit can compile
    without HDF5's C++ API. If `H5Cpp.h` is ever put behind a conditional, a
    build without HDF5 becomes conceivable and the FATAL_ERROR needs revisiting.
    """
    generic = (REPO / "Src" / "Generic.hpp").read_text()
    lines = generic.splitlines()

    # The file's own include guard wraps everything and must not count.
    guard = None
    for index, line in enumerate(lines):
        match = re.match(r"#\s*ifndef\s+(\w+)", line.strip())
        if match:
            rest = [l.strip() for l in lines[index + 1:index + 3]]
            if any(re.match(rf"#\s*define\s+{match.group(1)}\b", r) for r in rest):
                guard = index
            break

    depth, found = 0, None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if re.match(r"#\s*(if|ifdef|ifndef)\b", stripped):
            if index != guard:
                depth += 1
        elif re.match(r"#\s*endif\b", stripped):
            depth = max(0, depth - 1)
        elif re.match(r"#\s*include\b", stripped) and "H5Cpp.h" in stripped:
            found = (stripped, depth)
            break

    assert found is not None, "Src/Generic.hpp no longer includes H5Cpp.h"
    text, depth = found
    assert depth == 0, (
        f"{text!r} sits inside a conditional ({depth} levels deep, excluding the "
        "include guard); the FATAL_ERROR in CMakeLists.txt needs revisiting.")
