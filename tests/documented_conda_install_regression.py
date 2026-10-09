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


def _dependencies_raw():
    """The dependency specs in environment.yml, with their bounds.

    Parsed with a regex rather than a YAML library: PyYAML is importable in the
    development environment but is deliberately not a declared dependency, and
    adding one so a test can read a 13-line file is the wrong trade. The format
    is a flat list of scalars; a nested `pip:` sub-list would need revisiting.
    """
    body = ENVIRONMENT.read_text().split("dependencies:", 1)[1]
    specs = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped.startswith("- "):
            continue
        specs.append(stripped[2:].split("#", 1)[0].strip())
    return specs


def _dependencies():
    """Just the package names, without bounds or build strings."""
    return [re.split(r"[ >=<!~]", spec, maxsplit=1)[0] for spec in _dependencies_raw()]


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


#: A real lower bound: a comparison operator followed by a digit. Testing for
#: the mere presence of `<>=!~` is not enough -- `hdf5=*=nompi*` contains `=`
#: and `numpy=` contains `=`, and neither constrains anything. conda's own
#: MatchSpec agrees on the second but not the first: it reports version `'*'`
#: for the wildcard, which is indistinguishable from a bound if you only check
#: that a version component exists.
LOWER_BOUND = re.compile(r">=\s*\d")


def test_every_conda_dependency_carries_a_lower_bound():
    """An unconstrained dependency is a silent upgrade path.

    Ten of the thirteen entries once had no bound at all, while the page called
    the file reproducible. The bounds are floors, not pins -- the point is to
    state what KITE needs, not to freeze a solve.
    """
    unbounded = [spec for spec in _dependencies_raw()
                 if not LOWER_BOUND.search(spec)]
    assert not unbounded, (
        "these dependencies have no lower bound, so the documented install "
        f"silently accepts any release: {unbounded}")


def test_every_pip_dependency_carries_a_lower_bound():
    """pyproject.toml is the other half, and the only half three routes get.

    `pip install -e .` ends every documented install route, conda or not, so
    bounding environment.yml alone leaves the Linux, macOS and Docker routes
    with no Python floors at all.
    """
    body = (REPO / "pyproject.toml").read_text().split("dependencies = [", 1)[1]
    body = body.split("]", 1)[0]
    specs = [line.strip().strip('",') for line in body.splitlines()
             if line.strip().startswith('"')]
    assert specs, "pyproject.toml declares no dependencies; has the format changed?"
    unbounded = [s for s in specs if not LOWER_BOUND.search(s)]
    assert not unbounded, (
        f"these pip dependencies have no lower bound: {unbounded}")


def test_the_two_dependency_files_agree_on_python():
    """One file permitting 3.13 while the other forbids it is a real ambiguity."""
    pyproject = (REPO / "pyproject.toml").read_text()
    match = re.search(r'requires-python\s*=\s*"([^"]+)"', pyproject)
    assert match, "pyproject.toml no longer declares requires-python"
    conda = next(s for s in _dependencies_raw() if s.startswith("python"))
    conda_range = conda.split(None, 1)[1].strip()
    assert match.group(1) == conda_range, (
        f"pyproject.toml allows Python {match.group(1)!r} but environment.yml "
        f"allows {conda_range!r}; pip and conda would disagree about 3.13")


def test_the_page_does_not_claim_the_environment_is_reproducible():
    """`environment.yml` is a set of floors; it is not a lockfile.

    Calling it reproducible promised something it cannot deliver -- two solves
    weeks apart give different versions, as the bounds above allow by design.
    Seed-based reproducibility in docs/api/kite.md is a different claim and is
    untouched by this.
    """
    page = PAGE.read_text()
    assert "{ #conda_section }" in page, (
        "the conda section is gone, so this test is guarding nothing")

    heading = next(line for line in page.splitlines()
                   if "{ #conda_section }" in line)
    assert "reproducib" not in heading.lower(), (
        f"the conda section heading claims reproducibility again: {heading!r}. "
        "environment.yml carries lower bounds only; back the claim with a "
        "lockfile or do not make it.")

    # Scoped to the section, not the whole page: a truthful sentence saying the
    # file is NOT a lockfile is a documentation improvement, and a blanket ban
    # on the word would forbid it.
    section = page.split("{ #conda_section }", 1)[1].split("\n## ", 1)[0]
    claims = [line.strip() for line in section.splitlines()
              if "reproducib" in line.lower() and "not " not in line.lower()]
    assert not claims, (
        f"section 2.4 claims reproducibility: {claims}")


def test_the_hdf5_nompi_pin_survives_alongside_its_version_bound():
    """Adding a version bound must not drop the build-string pin.

    `hdf5 >=1.10 nompi*` is name + version range + build glob. Writing it as
    `hdf5 >=1.10` would resolve to an MPI build and drag in mpich/openmpi.
    """
    raw = ENVIRONMENT.read_text()
    match = re.search(r"-\s*hdf5\s+([^\n#]+)", raw)
    assert match, "hdf5 is no longer declared"
    spec = match.group(1).strip()
    assert "nompi" in spec, f"the nompi build pin is gone from {spec!r}"
    assert re.search(r"[<>=]", spec), f"hdf5 has no version bound in {spec!r}"


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
