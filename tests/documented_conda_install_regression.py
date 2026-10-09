"""The conda install page and environment.yml must keep saying the same thing.

docs/installation.md section 2.4 is the recommended install route. It names the
environment, promises what the file provides, and prints the exact commands.
None of that is checked anywhere, so a package dropped from `environment.yml`
or an environment renamed leaves the page quietly wrong -- and the reader only
finds out partway through an install.

This does not build anything. Whether the route actually produces working
binaries is covered by CI's `build-and-test` job, which runs the same file on
ubuntu-latest and macos-latest; verified by hand on macOS as well, where it
configures against conda's HDF5 2.1.0, FFTW3 and OpenMP 5.1, builds both
binaries, and computes a graphene DOS end to end.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ENVIRONMENT = REPO / "environment.yml"
PAGE = REPO / "docs" / "installation.md"


def _dependencies():
    """The dependency names in environment.yml, without their version pins."""
    text = ENVIRONMENT.read_text()
    body = text.split("dependencies:", 1)[1]
    names = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped.startswith("- "):
            continue
        names.append(re.split(r"[ >=<=!~]", stripped[2:], maxsplit=1)[0])
    return names


def _env_name():
    match = re.search(r"^name:\s*(\S+)", ENVIRONMENT.read_text(), re.M)
    assert match, "environment.yml has no name:"
    return match.group(1)


def test_the_page_activates_the_environment_the_file_declares():
    name = _env_name()
    page = PAGE.read_text()
    assert f"conda activate {name}" in page, (
        f"environment.yml declares {name!r} but the page does not activate it")


@pytest.mark.parametrize("promise,packages", [
    ("a matched C/C++ compiler toolchain", ("c-compiler", "cxx-compiler")),
    ("HDF5", ("hdf5",)),
    ("FFTW3", ("fftw",)),
])
def test_the_file_provides_what_the_page_promises(promise, packages):
    """Both sides: the page must make the promise, and the file must keep it.

    Asserting only against environment.yml would make this a snapshot of that
    file -- it would pass unchanged if section 2.4 were deleted outright, and
    it would "detect" nothing beyond its own mutation.
    """
    page = PAGE.read_text()
    assert promise in page, (
        f"the page no longer promises {promise!r}; this test is pinning a claim "
        "that is not made any more, so update or remove it")
    declared = _dependencies()
    for package in packages:
        assert package in declared, (
            f"the page promises {promise!r} but environment.yml does not list "
            f"{package!r}; the recommended install route would not deliver it")


def test_openmp_is_in_the_environment_because_the_build_now_requires_it():
    """CMake fails hard without OpenMP, so the recommended route must supply it.

    The second oracle here is CMakeLists.txt rather than the page: the guard is
    what creates the requirement. See tests/openmp_required_regression.py.
    """
    cmake = (REPO / "CMakeLists.txt").read_text()
    assert "OPENMP_FOUND" in cmake and "FATAL_ERROR" in cmake, (
        "CMake no longer fails without OpenMP, so this test's premise is gone")
    assert "llvm-openmp" in _dependencies(), (
        "CMake fails hard without OpenMP but environment.yml no longer supplies "
        "it; the documented install would produce a tree that cannot configure")


def test_hdf5_is_pinned_to_a_build_without_mpi():
    """`hdf5=*=nompi*` is a build-string pin, not a version pin.

    The MPI builds pull in a parallel HDF5 that KITE does not use and that
    drags mpich/openmpi into the environment. Recorded because the pin looks
    like a typo and is easy to 'tidy' away.
    """
    raw = ENVIRONMENT.read_text()
    assert re.search(r"-\s*hdf5=\*=nompi\*", raw), (
        "the nompi build pin on hdf5 is gone; the environment may now resolve "
        "to a parallel HDF5 build")


def test_the_micromamba_route_initialises_the_shell_before_activating():
    """`micromamba activate` without the shell hook fails outright.

    Measured: `critical libmamba Shell not initialized ... micromamba is running
    as a subprocess and can't modify the parent shell`. The page handles this
    for conda (`conda init`) and once omitted it for micromamba, which stopped
    a MacPorts reader at the third command of the recommended route.
    """
    page = PAGE.read_text()
    if "micromamba activate" not in page:
        pytest.skip("the page no longer documents a micromamba route")
    hook = page.index("micromamba shell hook")
    activate = page.index("micromamba activate")
    assert hook < activate, (
        "`micromamba activate` is documented before the shell hook that makes "
        "it work; following the page in order fails")
