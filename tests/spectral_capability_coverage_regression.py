"""Each spectral capability names the test that protects it, and that test passes.

RF-001 to RF-006 are the capabilities this cycle does not change: the density of
states, the local density of states, the real-space LDoS map, ARPES, the
momentum-space map, and the band structure computed without the KPM machinery.
"Protected by the legacy suite" is only true if someone can say *which* directory
protects *which* capability -- otherwise a directory can quietly stop covering
anything and the claim survives, which is how a suite rots while still reporting
success.

So this file is the mapping, written down and executed. For the capabilities the
legacy suite does cover it runs those directories and checks their verdict; for
the two it never covered -- the LDoS map and the momentum map, both added this
cycle -- it names the modern test instead and checks that it is still there.

RF-001 says "in 2D and 3D", and exactly one directory carries the 3D half. That is
asserted against the stored configuration rather than the filename, because
`cb` meaning "cubic" is a convention, not a check.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TESTS = REPO / "tests"
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

#: capability -> (legacy directories, modern test files)
COVERAGE = {
    "RF-001": (["test001_KITEx_sq", "test002_KITEx_gr", "test003_KITEx_ha",
                "test004_KITEx_cb", "test017_KITEt_gr_opt",
                "test018_KITEt_gr_doslim"],
               []),
    "RF-002": (["test008_KITEx_sq_ldos", "test011_KITEx_sq_defec",
                "test014_KITEx_gr_adat"],
               ["operator_ldos_regression.py",
                "operator_ldos_analytic_stochastic_regression.py"]),
    "RF-003": ([], ["ldos_map_reader_regression.py"]),
    "RF-004": (["test010_KITEx_sq_arpes", "test012_KITEx_gr_arpes",
                "test013_KITEx_km_arpes"],
               ["arpes_exclusive_regression.py"]),
    "RF-005": ([], ["spectral_map_reader_regression.py"]),
    "RF-006": ([], ["hamiltonian_k_convention_regression.py"]),
}

#: The one directory carrying RF-001's three-dimensional half.
THREE_DIMENSIONAL = "test004_KITEx_cb"

LEGACY = sorted({d for dirs, _ in COVERAGE.values() for d in dirs})


@pytest.mark.parametrize("capability", sorted(COVERAGE))
def test_every_capability_names_something(capability):
    legacy, modern = COVERAGE[capability]
    assert legacy or modern, f"{capability} names no test at all"
    for name in legacy:
        assert (TESTS / name / "test.py").exists(), (
            f"{capability} names {name}, which has no test.py")
    for name in modern:
        assert (TESTS / name).exists(), (
            f"{capability} names {name}, which is not in tests/")


def test_the_density_of_states_is_covered_in_three_dimensions():
    """RF-001 says 2D *and* 3D, and only one directory carries the 3D half."""
    h5py = pytest.importorskip("h5py")
    path = TESTS / THREE_DIMENSIONAL / "configORIG.h5"
    assert path.exists(), f"{THREE_DIMENSIONAL} has no stored configuration"
    with h5py.File(path, "r") as f:
        dim = int(f["DIM"][()])
        assert "/Calculation/dos" in f, (
            f"{THREE_DIMENSIONAL} no longer requests a density of states, so "
            "nothing covers RF-001 in three dimensions")
    assert dim == 3, (
        f"{THREE_DIMENSIONAL} is {dim}-dimensional. RF-001 requires a 3D case "
        "and this is the only directory that provided one")


@pytest.mark.skipif(not (KITEX.exists() and KITE_TOOLS.exists()),
                    reason="needs a built KITEx and KITE-tools")
@pytest.mark.parametrize("directory", LEGACY)
def test_the_named_directory_passes(directory):
    """The quick path: restore the stored input, run, read the verdict.

    `config.h5` is gitignored while `configORIG.h5` and `configREF.h5` are
    tracked, so running in place leaves the working tree clean.
    """
    d = TESTS / directory
    shutil.copy(d / "configORIG.h5", d / "config.h5")

    env = dict(os.environ, KITEX=str(KITEX), KITE_TOOLS=str(KITE_TOOLS))
    proc = subprocess.run([sys.executable, "test.py"], cwd=d, env=env,
                          capture_output=True, text=True)
    out = proc.stdout

    assert proc.returncode == 0, f"{directory}/test.py exited {proc.returncode}"
    assert "ERROR" not in out, (
        f"{directory} could not run KITEx: {out.strip()!r}")
    assert "Problem" not in out, (
        f"{directory} disagrees with its stored reference: {out.strip()!r}")
    assert "OK" in out, (
        f"{directory} printed no verdict at all: {out.strip()!r}. A test that "
        "says nothing is indistinguishable from one that passed")
