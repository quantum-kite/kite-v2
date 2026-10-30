"""`custom_two`'s moment count must be checked against the block size too.

Same defect as `Gamma2D`/`Gamma3D`, in `SimulationRankTwo.cpp` and
`SimulationRankTwoLocal.cpp`: the guard checked `% 2` while the accumulation
steps by `MEMORY` (16) and runs a whole block whatever is left, so
`ind = m + i + M*(n + j)` walks past `gamma`.

When this was written nothing on the Python side validated it, so an ordinary
`custom.Vertex(200, ...)` reached the guard through the normal API; the exporter
now refuses it too, so the file here is patched after export, which is what a
hand-authored one amounts to. Demonstrated by instrumenting the write -- 1672
writes past the end of a 40000-element array, the worst at index 41607; at 192
moments, none. It does not crash, so the moments are simply wrong.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

MEMORY = 16
VALID = 192
INVALID = 200        # even, but 200 % 16 == 8

_BUILDER = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

M = int(sys.argv[2])
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['up', [0, 0], 0.0], ['dn', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'up', 'up', -1.0), ([0, 1], 'up', 'up', -1.0),
                 ([1, 0], 'dn', 'dn', -1.0), ([0, 1], 'dn', 'dn', -1.0),
                 ([0, 0], 'up', 'dn', 0.2))
cfg = kite.Configuration(divisions=[1, 1], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-6, 6])
calc = kite.Calculation(cfg)
calc.add_orbital_index('up', 0)
calc.add_orbital_index('dn', 1)
calc.add_orbital_coupling('up', 'up', 0.5, 'l0')
calc.add_orbital_coupling('dn', 'dn', -0.5, 'l0')
calc.custom_two(stream_=[custom.Vertex(M, [[0.5, 'vx.l0'], [0.5, 'l0.vx']]),
                         custom.Vertex(M, [[1.0, 'vy']])],
                num_random_=1, num_disorder_=1, num_points_=100,
                temperature_=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _run(tmp, moments):
    """Export at a valid count, then patch it -- the exporter refuses the bad
    one now, and this guard is the backstop for files it never saw."""
    h5py = pytest.importorskip("h5py")
    path = tmp / f"c{moments}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path), str(VALID)],
                   cwd=tmp, check=True, capture_output=True, text=True)
    if moments != VALID:
        with h5py.File(path, "r+") as f:
            for vertex in ("Vertex0", "Vertex1"):
                f[f"/Calculation/CustomTwo/{vertex}/NumMoments"][...] = moments
    return subprocess.run([str(KITEX), str(path)], capture_output=True, text=True)


@pytest.fixture(scope="module")
def rejected(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("ct_bad"), INVALID)


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("ct_ok"), VALID)


def test_the_file_reaches_the_guard_at_all(rejected):
    """If the run never got as far as the guard, the rest proves nothing."""
    assert "Traceback" not in (rejected.stdout + rejected.stderr)


def test_an_even_count_that_is_not_a_multiple_of_the_block_is_rejected(rejected):
    assert rejected.returncode == 1, (
        f"{INVALID} moments is even, so the old `% 2` guard passed it and the "
        "accumulation wrote past gamma")


def test_the_message_names_the_block_size(rejected):
    lines = [l for l in (rejected.stdout + rejected.stderr).splitlines()
             if "number of moments" in l.lower()]
    assert lines, "the guard printed no message"
    assert str(MEMORY) in lines[0], (
        f"the message does not say what to change {INVALID} to: {lines[0]!r}")


def test_a_multiple_of_the_block_size_still_runs(accepted):
    assert accepted.returncode == 0, (
        f"{VALID} is a multiple of {MEMORY} and must be accepted")
    assert "Done" in accepted.stdout
