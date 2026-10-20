"""A file without /Hamiltonian/CustomLocalPotential must mean "no custom potential".

KITEx reads that dataset into `check_local_potential` inside a `try` whose
`catch (...)` is empty (Src/Hamiltonian/Hamiltonian.cpp). The dataset is
optional -- the exporter only began writing it in April 2026 -- so on any older
file the read throws and the flag keeps whatever was on the stack. When that
happened to be non-zero, `Hamiltonian.cpp` added the demo potential from
`Src/custom.cpp`, which saturates at 2.0, to every Anderson on-site entry. At
this fixture's energy scale that is +0.5 rescaled: the spectrum leaves [-1, 1]
and the Chebyshev recursion diverges.

It fails open rather than closed, and silently: no error, no warning, a result
of order 1e32 where a conductivity of order 10 belongs.

All eighteen legacy fixtures predate the dataset, and this one is the only place
it shows, because the spurious term is only read when an orbital carries an
Anderson model. The flag is now initialised to 0 in the constructor.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
#: Carries Anderson disorder and predates the dataset, so it exercises the path.
FIXTURE = REPO / "tests" / "test007_KITEx_gr_ss" / "configORIG.h5"
RESULT = "/Calculation/singleshot_conductivity_dc/SingleShot"

needs_build = pytest.mark.skipif(
    not (KITEX.exists() and FIXTURE.exists()),
    reason="needs a built KITEx and the test007 fixture")


def _run(where, add_dataset):
    where.mkdir(parents=True, exist_ok=True)
    local = where / "config.h5"
    shutil.copy(FIXTURE, local)
    with h5py.File(local, "r+") as handle:
        present = "Hamiltonian/CustomLocalPotential" in handle
        if add_dataset and not present:
            handle.create_dataset("Hamiltonian/CustomLocalPotential",
                                  data=np.uint32(0), dtype="u4")
        elif not add_dataset:
            assert not present, "the fixture now carries the dataset; pick another"
    done = subprocess.run([str(KITEX), local.name], cwd=str(where),
                          capture_output=True, text=True)
    assert done.returncode == 0, f"KITEx failed:\n{done.stderr[-400:]}"
    with h5py.File(local, "r") as handle:
        return np.array(handle[RESULT])[:, 3]


@needs_build
def test_a_file_without_the_dataset_does_not_diverge(tmp_path):
    sigma = _run(tmp_path, add_dataset=False)
    assert np.isfinite(sigma).all(), "the conductivity is not finite"
    assert np.abs(sigma).max() < 1e3, (
        f"max |sigma| is {np.abs(sigma).max():.3e}; a missing optional dataset is "
        "being read as a request for the custom potential, so the rescaled "
        "spectrum leaves [-1, 1] and the Chebyshev recursion diverges")


@needs_build
def test_present_and_absent_give_the_same_answer(tmp_path):
    """Absent must mean zero, not merely 'something finite'."""
    without = _run(tmp_path / "a", add_dataset=False)
    with_zero = _run(tmp_path / "b", add_dataset=True)
    # This fixture is PRECISION 0, i.e. float32, and the moments are summed
    # across threads in whatever order they finish, so agreement floors at
    # about 1e-7 relative. The defect this guards against is a factor of 1e31.
    assert without == pytest.approx(with_zero, rel=1e-5), (
        "a file lacking the dataset gives a different answer from one carrying "
        "an explicit zero, so the default is not 'no custom potential'")
