"""Configuring without OpenMP must fail and say what to do; with it, succeed.

`#include <omp.h>` is unconditional in `Src/Generic.hpp` and in twelve
translation units under `tools/Src/`, so there is no serial build to fall back
to. CMake's `find_package(OpenMP)` without `REQUIRED` used to let configuration
succeed anyway, and the missing dependency surfaced as a wall of C++ errors in
every file instead of one line at configure time.

Both directions are covered deliberately. A guard that never fires is useless;
a guard that *always* fires is worse, because it kills every build on every
platform, and a suite that only exercises the failing path cannot tell the two
apart. Whether OpenMP is actually available is established with a throwaway
project rather than with KITE's own CMakeLists, so a broken guard cannot
masquerade as a runner without OpenMP.

OpenMP is suppressed with `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON`, a standard
CMake variable, rather than by removing a compiler -- the point is to exercise
the project's own branch, not the toolchain's.
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
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "PROBE_OPENMP_YES")
endif()
"""


def _configure(source, build, *extra):
    return subprocess.run(["cmake", "-S", str(source), "-B", str(build), *extra],
                          capture_output=True, text=True)


@pytest.fixture(scope="module")
def without_openmp(tmp_path_factory):
    """One configure, shared: it costs seconds and both tests read the same run."""
    build = tmp_path_factory.mktemp("no-omp")
    return _configure(REPO, build, "-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON")


@pytest.fixture(scope="module")
def openmp_is_available(tmp_path_factory):
    """Does this toolchain have OpenMP at all?

    Answered with a minimal project of our own, NOT with KITE's CMakeLists. If
    we asked KITE's, a guard broken so that it always fires would look exactly
    like a machine without OpenMP, and the positive test below would skip
    itself into uselessness.
    """
    root = tmp_path_factory.mktemp("probe")
    (root / "CMakeLists.txt").write_text(PROBE)
    proc = _configure(root, root / "build")
    return "PROBE_OPENMP_YES" in proc.stdout


@needs_cmake
def test_configure_fails_without_openmp_and_says_what_to_do(without_openmp):
    """Both halves together: it must fail, and the failure must be actionable.

    Checked as one test because a non-zero exit on its own proves little -- a
    missing FFTW would also produce it.
    """
    assert without_openmp.returncode != 0, (
        "configuration succeeded without OpenMP; the build would then fail in "
        "every translation unit instead:\n" + without_openmp.stdout)
    message = without_openmp.stdout + without_openmp.stderr
    assert "OpenMP" in message
    for expected in ("omp.h", "CMAKE_CXX_COMPILER", "environment.yml",
                     "docs/installation.md"):
        assert expected in message, (
            f"the failure does not mention {expected!r}, so a user cannot act on it:\n"
            + message)


@needs_cmake
def test_configure_succeeds_when_openmp_is_present(openmp_is_available, tmp_path):
    """The other direction: the guard must not fire when it should not.

    Without this, a guard rewritten to fire unconditionally -- killing every
    build everywhere -- leaves the suite entirely green.
    """
    if not openmp_is_available:
        pytest.skip("this toolchain has no OpenMP, so the positive path cannot be tested")
    proc = _configure(REPO, tmp_path / "with-omp")
    assert proc.returncode == 0, (
        "OpenMP is available on this toolchain but configuration still failed:\n"
        + proc.stdout + proc.stderr)
    flags = (tmp_path / "with-omp" / "CMakeFiles" / "KITEx.dir" / "flags.make").read_text()
    assert "-fopenmp" in flags or "-openmp" in flags or "/openmp" in flags, (
        "configuration succeeded but no OpenMP flag reached the compile line; the "
        f"set(...) lines after the guard may have been dropped:\n{flags[:400]}")
    assert "-O3" in flags, "the -O3 that the OpenMP branch used to add is gone"


def test_generic_hpp_includes_omp_h_unconditionally():
    """The guard's premise, read from the source.

    Failing hard is only justified while no file can compile without omp.h. If
    the include is ever put behind `#ifdef _OPENMP`, a serial build becomes
    possible and the FATAL_ERROR should be reconsidered -- so check that the
    include is reached unconditionally, not merely that it is present.
    """
    lines = (REPO / "Src" / "Generic.hpp").read_text().splitlines()

    # The file's own include guard is `#ifndef NAME` followed by `#define NAME`.
    # It wraps everything and means nothing here, so it must not count as a
    # conditional -- otherwise this check fires on every well-formed header.
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
        elif re.match(r"#\s*include\b", stripped) and "omp.h" in stripped:
            found = (stripped, depth)
            break

    assert found is not None, "Src/Generic.hpp no longer includes <omp.h>"
    text, depth = found
    assert depth == 0, (
        f"{text!r} sits inside a conditional ({depth} levels deep, excluding the "
        "include guard). A serial build may now be possible, so the FATAL_ERROR "
        "in CMakeLists.txt needs revisiting.")
