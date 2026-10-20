"""`-X` means only that quantity, for ARPES as for everything else.

KITE-tools decides what to compute in two steps: count how many `-X` flags were
given, then either honour the one exclusive quantity or compute everything.
`ARPES_Exclusive` was set by the parser but never counted, so `--ARPES -X`
produced a count of zero, fell into the "compute everything" branch, and ran
every other quantity in the file as well. A file carrying a single calculation
could never show it -- there was nothing else to wrongly compute.

The fixture therefore carries two calculations, ARPES and DOS.
"""

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
H5 = FIXTURES / "square_arpes_and_dos_check-output.h5"
KITE_TOOLS = REPO / "build" / "KITE-tools"

needs_fixture = pytest.mark.skipif(
    not (H5.exists() and KITE_TOOLS.exists()),
    reason="needs tests/fixtures/square_arpes_and_dos_check-output.h5 and a built "
           "KITE-tools; run the fixture script and KITEx first")


def _run(tmp_path, *flags):
    import shutil
    local = tmp_path / H5.name
    shutil.copy(H5, local)
    subprocess.run([str(KITE_TOOLS), local.name, *flags],
                   cwd=str(tmp_path), capture_output=True, text=True, check=True)
    return {p.name for p in tmp_path.glob("*.dat")}


@needs_fixture
def test_arpes_exclusive_skips_the_density_of_states(tmp_path):
    produced = _run(tmp_path, "--ARPES", "-X")
    assert "arpes.dat" in produced, "ARPES was requested but not produced"
    assert "dos.dat" not in produced, (
        "--ARPES -X also computed the density of states; the exclusive flag was "
        f"not honoured (produced {sorted(produced)})")


@needs_fixture
def test_dos_exclusive_skips_arpes(tmp_path):
    """The other direction, to show the fixture can tell the two apart."""
    produced = _run(tmp_path, "--DOS", "-X")
    assert "dos.dat" in produced
    assert "arpes.dat" not in produced


@needs_fixture
def test_without_any_exclusive_flag_both_are_computed(tmp_path):
    produced = _run(tmp_path)
    assert {"arpes.dat", "dos.dat"} <= produced, (
        f"a bare run should compute everything in the file (produced {sorted(produced)})")
