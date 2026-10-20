"""KITE-tools reports the single-shot conductivity at the energies asked for.

Single-shot is the one calculation KITEx finishes by itself: it solves the
linear system at each Fermi energy, so there is no Chebyshev series left to
reconstruct and nothing for KITE-tools to choose. `--SingleShot` therefore only
has to report what is stored, at the right energies.

"The right energies" is the whole content of this file. KITE writes energies
rescaled to the Chebyshev interval, `(E - shift)/scale`, so reporting them
again means undoing both. The fixture uses an asymmetric spectrum range on
purpose: with a symmetric one the shift is zero and dropping it is invisible.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
H5 = FIXTURES / "graphene_singleshot_check-output.h5"
KITE_TOOLS = REPO / "build" / "KITE-tools"
HELPER = REPO / "tools" / "process_single_shot.py"

sys.path.insert(0, str(FIXTURES))
from graphene_singleshot_check import ENERGIES  # noqa: E402

needs_fixture = pytest.mark.skipif(
    not (H5.exists() and KITE_TOOLS.exists()),
    reason="needs tests/fixtures/graphene_singleshot_check-output.h5 and a built "
           "KITE-tools; run the fixture script and KITEx first")


def _kite_tools(tmp_path):
    """Run `--SingleShot` on a copy, since the .dat lands beside the input."""
    import shutil
    local = tmp_path / H5.name
    shutil.copy(H5, local)
    subprocess.run([str(KITE_TOOLS), local.name, "--SingleShot", "-X"],
                   cwd=str(tmp_path), capture_output=True, text=True, check=True)
    return np.loadtxt(tmp_path / "singleshot_condDC.dat")


@needs_fixture
def test_the_energies_come_back_as_they_were_asked_for(tmp_path):
    data = _kite_tools(tmp_path)
    assert data[:, 0] == pytest.approx(ENERGIES, abs=1e-6), (
        "the reported Fermi energies are not the ones the calculation was set "
        "up with; the rescaling to the Chebyshev interval has not been undone")


@needs_fixture
def test_the_conductivity_is_the_number_kitex_stored(tmp_path):
    """Nothing is recomputed, so only the .dat formatting may differ.

    KITE-tools writes every quantity at the stream default of six significant
    figures, so the tolerance here is that of the file format, not of a
    calculation -- there is no calculation.
    """
    import h5py
    data = _kite_tools(tmp_path)
    with h5py.File(H5, "r") as handle:
        stored = np.asarray(
            handle["/Calculation/singleshot_conductivity_dc/SingleShot"])[:, 3]
    assert data[:, 1] == pytest.approx(stored, rel=1e-5), (
        "KITE-tools is not reporting the conductivity KITEx computed")


@needs_fixture
def test_the_python_helper_agrees_on_the_conductivity_but_not_the_energy(tmp_path):
    """Pins the one difference between the two routes, so it is not a surprise.

    `tools/process_single_shot.py` reads the stored table directly, and that
    table's energy column carries the energy scale but not the energy shift.
    The conductivities are the same numbers; the energies differ by the shift.
    Pinned here so that if the stored column is ever corrected, the test fails
    and says so rather than the two routes drifting apart unnoticed.
    """
    import shutil
    import h5py
    local = tmp_path / H5.name
    shutil.copy(H5, local)
    subprocess.run([sys.executable, str(HELPER), local.name],
                   cwd=str(tmp_path), capture_output=True, text=True, check=True)
    helper = np.loadtxt(tmp_path / (local.stem + ".dat"))
    tools = _kite_tools(tmp_path)

    with h5py.File(H5, "r") as handle:
        shift = float(np.asarray(handle["EnergyShift"]).ravel()[0])
    assert shift != 0, "the fixture must have an asymmetric range for this to bite"

    assert helper[:, 2] == pytest.approx(tools[:, 1], rel=1e-5)
    assert helper[:, 0] == pytest.approx(tools[:, 0] - shift, abs=1e-6), (
        "the helper's energies no longer differ from KITE-tools' by exactly the "
        "energy shift; one of the two routes has changed convention")
