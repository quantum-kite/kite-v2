"""The paper scripts must still build an HDF5 that KITEx would accept.

`examples/paper/` reproduces the figures of the KITE paper. The scripts were
written against the sizes used there -- tens of thousands of unit cells -- so
nobody re-runs them casually, and a change to the Python exporter can break
them without anything noticing.

This checks the export only: each producer runs to `kite.config_system` and the
resulting file is opened and required to carry both groups KITEx reads. KITEx
itself is never invoked -- the physics of these figures is not under test here,
the file format is.

Sizes are cut where the script allows it, because the point is the file rather
than the run. Three scripts set their own floor and are given their published
size instead: a magnetic field has a minimum representable value set by the
system size, a wave packet needs a k-grid fine enough to contain a state at the
requested energy, and a vacancy pinned at site (32, 32) needs a lattice that
reaches it. Those are properties of the examples, not limits of this test.
"""

import contextlib
import io
import os
import runpy
import shutil
import sys
from pathlib import Path

import pytest

h5py = pytest.importorskip("h5py")

REPO = Path(__file__).resolve().parents[1]
PAPER = REPO / "examples" / "paper"

#: Arguments per producer. An empty list means the script's own defaults, which
#: are the published values -- used where a smaller system cannot be exported.
ARGV = {
    "anthracene.py": ["16"],
    "benzene.py": ["16"],
    "naphthalene.py": ["16"],
    "oligoacenes.py": ["16"],
    "pentacene.py": ["16"],
    "tetracene.py": ["16"],
    "qahe_disorder.py": ["16", "1", "1"],
    "dos_gr_nn.py": ["16", "16", "1", "1", "32"],
    "dos_gr_nnn.py": ["16", "16", "1", "1", "32"],
    "config.py": ["16", "32", "0.5"] + ["1"] * 8 + ["0", "0.0"],
    # The vacancy sits at (32, 32), so the lattice has to reach it.
    "config_ldos.py": ["64", "32"],
    # Below the published size the minimum representable field exceeds the
    # requested one and KITE declines to write the calculation.
    "phosphorene_magnetic.py": [],
}

#: These read an HDF5 or a .dat and plot it. They are checked by source rather
#: than by running: without the figure data they have nothing to open, so
#: running them would only prove that the input is absent.
CONSUMERS = {
    "dataplots2.py", "plot_benchmarks.py", "plot_dos.py",
    "plot_dos_jackson_splitted.py", "plot_spin.py", "process_arpes.py",
    "process_arpes1.py", "process_new.py",
}

#: Scans a k-grid for a state at the requested energy and finds none, so the
#: spinor comes out empty and the export stops. The window is set by
#: `k = energy / vf` with the energy in eV and vf in m/s, which is not a
#: wavevector; the k-space window is then ~1e-7 of the zone at every physical
#: energy. The line predates this work (079c510, 2022). Left as found: the
#: paper folder is out of scope, and the fix is a physics decision, not a
#: mechanical one.
BLOCKED = {
    "gaussian_wavepacket_only_anderson.py",
    "gaussian_wavepacket_magnetic_impurities.py",
}

#: Builds its lattice with pybinding, which KITE no longer depends on.
NEEDS_PYBINDING = {"tblg_DOS.py"}


def _scripts():
    return sorted(PAPER.rglob("*.py"), key=lambda p: p.name)


def _run(script, work):
    """Run `script` inside `work` with its arguments, muted."""
    for extra in script.parent.glob("*"):
        if extra.is_file() and not extra.name.endswith(".h5"):
            shutil.copy(extra, work / extra.name)
    cwd, argv = os.getcwd(), sys.argv[:]
    try:
        os.chdir(work)
        sys.argv = [script.name] + ARGV.get(script.name, [])
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            try:
                runpy.run_path(script.name, run_name="__main__")
            except SystemExit:
                pass
    finally:
        os.chdir(cwd)
        sys.argv = argv


@pytest.mark.parametrize("script", _scripts(), ids=lambda p: p.name)
def test_every_paper_script_exports_or_is_accounted_for(script, tmp_path, monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    name = script.name

    if name in NEEDS_PYBINDING:
        pytest.skip("builds its lattice with pybinding")
    if name in BLOCKED:
        pytest.xfail("empty spinor from the k-window; see BLOCKED in this file")

    if name in CONSUMERS:
        assert "config_system" not in script.read_text(), (
            f"{name} is catalogued as a consumer but calls kite.config_system; "
            "it exports, so it belongs in ARGV and must be run")
        return

    _run(script, tmp_path)
    written = sorted(tmp_path.glob("*.h5"))

    assert written, (
        f"examples/paper/{script.relative_to(PAPER)} ran but exported no HDF5; "
        "it is neither a consumer nor a known-blocked script")
    for path in written:
        with h5py.File(path, "r") as handle:
            missing = [g for g in ("Hamiltonian", "Calculation") if g not in handle]
        assert not missing, (
            f"{name} wrote {path.name} without {missing}; KITEx cannot read it. "
            "A partial file like this means the script stopped after the lattice "
            "but before the calculation was added")
