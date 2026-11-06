"""A hopping address that overflows the orbital field must stop the run.

The destination of a hopping is one integer: `orbital*3**D + sum_d (n_d + 1)*3**d`.
Each cell component gets one base-3 digit, so it can only be -1, 0 or +1. A
longer bond carries into the digit above, and where that digit is the orbital
field the decoded orbital runs past `Orb` -- `HamiltonianRegular.cpp` then
indexes `rOrb.col()` off the end of the table.

The rule originally said to warn and proceed. That was amended to abort
(PM, 2026-10-30) once proceeding was shown to be non-deterministic: with the
seeds fixed, six runs of the same file gave exit codes 1 1 1 1 1 0, because the
garbage read back varies.

What this canNOT catch, and the test says so rather than implying otherwise: a
carry out of a non-final direction produces a legitimate address for a
different cell. On this lattice `[+2,0]` and `[-1,+1]` both encode as 15, and
KITEx builds the second one -- correctly, as far as it can tell. Python refuses
both at export; this guard is only the backstop for hand-written files.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

# Two sublattices, four first-neighbour bonds. Seeded, so two runs of the same
# file are comparable moment by moment.
extra = sys.argv[2]
bonds = {'base':    [-1, 0],
         'collide': [-1, 1]}[extra]
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0], ['B', [0.5, 0.5], 0.0])
lat.add_hoppings(([0, 0] if extra == 'base' else bonds, 'A', 'B', -1.0),
                 ([-1, 0], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0),
                 ([-1, -1], 'A', 'B', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[16, 16],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[-5, 5],
                         seed_h=1234, seed_v=1234)
calc = kite.Calculation(cfg)
calc.dos(num_points=512, num_moments=64, num_random=1, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _export(tmp, name, flavour):
    path = tmp / f"{name}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path), flavour],
                   cwd=tmp, check=True, capture_output=True, text=True)
    return path


def _patch(path, a_to_b, b_to_a):
    """Overwrite the first A->B bond and its Hermitian partner."""
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "r+") as f:
        d = f["Hamiltonian/d"]
        a = np.asarray(d)
        a[0, 0] = a_to_b
        a[1, 3] = b_to_a
        d[...] = a


def _moments(path):
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "r") as f:
        g = f["/Calculation/dos"]
        return np.asarray(g["MU"][()]).ravel() if "MU" in g else None


def test_an_address_past_the_last_orbital_aborts(tmp_path):
    """[0,+2] carries into the orbital field: orbital 2 on a 2-orbital lattice."""
    path = _export(tmp_path, "orbflow", "base")
    _patch(path, 1 * 9 + (0 + 1) + 3 * (2 + 1), 0 * 9 + (0 + 1) + 3 * (-2 + 1))

    codes = [subprocess.run([str(KITEX), str(path)],
                            capture_output=True).returncode for _ in range(5)]
    assert set(codes) == {1}, (
        f"exit codes {codes}. Proceeding reads rOrb.col() past the table, and "
        "what comes back varies between runs")

    proc = subprocess.run([str(KITEX), str(path)], capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    assert "orbital 2" in out and "only 2" in out, (
        f"the message must name the hopping and the orbital: {out.strip()!r}")
    assert _moments(path) is None, "it must stop before writing any moments"


def test_a_valid_file_is_untouched(tmp_path):
    path = _export(tmp_path, "base", "base")
    proc = subprocess.run([str(KITEX), str(path)], capture_output=True, text=True)
    assert proc.returncode == 0
    assert "decodes to orbital" not in proc.stdout
    assert _moments(path) is not None


def test_a_carry_into_another_cell_digit_is_undetectable(tmp_path):
    """Not a gap to be fixed later: [+2,0] and [-1,+1] are the same integer.

    Asserted so that nobody later "fixes" this by adding a check that cannot
    exist, and so the limit of the guard above is on the record.
    """
    patched = _export(tmp_path, "patched", "base")
    _patch(patched, 1 * 9 + (2 + 1) + 3 * (0 + 1), 0 * 9 + (-2 + 1) + 3 * (0 + 1))
    collide = _export(tmp_path, "collide", "collide")

    for p in (patched, collide):
        assert subprocess.run([str(KITEX), str(p)],
                              capture_output=True).returncode == 0

    asked_for_two_cells = _moments(patched)
    one_cell_diagonal = _moments(collide)
    assert asked_for_two_cells is not None
    assert np.abs(asked_for_two_cells - one_cell_diagonal).max() == 0.0, (
        "a [+2,0] bond must decode to exactly the [-1,+1] model -- both encode "
        "as 15, so KITEx cannot tell them apart and builds the second")
