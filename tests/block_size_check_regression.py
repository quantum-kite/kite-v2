"""The moment count must be checked against the block size it is actually used with.

`Gamma2D` and `Gamma3D` accumulate in blocks of `MEMORY` (16): the outer loops
step `n += MEMORY` and the inner ones run the whole block regardless of what is
left. So a count that is not a multiple of `MEMORY` writes past `gamma`. The
guard checked `% 2`, which lets that through.

Demonstrated before the fix by instrumenting the write in `Gamma2D.cpp`: at 200
moments -- even, so the old guard passed it -- the flat index reached 41415 on a
40000-element array. At 192 there were no out-of-range writes at all. It does not
segfault on every machine, because the overrun is often absorbed by the
allocation; the result is silently wrong instead, which is worse.

The Python exporter already rejects this with a better message, so reaching the
C++ guard means a hand-authored HDF5 -- which is exactly what the guard is for.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

MEMORY = 16
VALID = 192          # 12 * 16
INVALID = 200        # even, but 200 % 16 == 8

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', -1.0), ([0, 1], 'A', 'A', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[-4.5, 4.5])
calc = kite.Calculation(cfg)
calc.conductivity_dc(direction='xy', num_points=100, num_moments=192,
                     num_random=1, num_disorder=1, temperature=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _run(tmp, moments):
    h5py = pytest.importorskip("h5py")
    path = tmp / f"m{moments}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path)], cwd=tmp,
                   check=True, capture_output=True, text=True)
    with h5py.File(path, "r+") as f:
        f["/Calculation/conductivity_dc/NumMoments"][...] = moments
    return subprocess.run([str(KITEX), str(path)], capture_output=True, text=True)


@pytest.fixture(scope="module")
def rejected(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("blocksize_bad"), INVALID)


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("blocksize_ok"), VALID)


def test_an_even_count_that_is_not_a_multiple_of_the_block_is_rejected(rejected):
    assert rejected.returncode == 1, (
        f"{INVALID} moments is even, so the old `% 2` guard let it through and "
        "the accumulation overran gamma. It must be refused")


def test_the_message_names_the_block_size_not_evenness(rejected):
    out = rejected.stdout + rejected.stderr
    # The licence banner contains the word "even", so look at the guard's own
    # line rather than the whole output.
    lines = [l for l in out.splitlines() if "number of moments" in l.lower()]
    assert lines, f"the guard printed no message at all: {out.strip()!r}"
    line = lines[0]
    assert str(MEMORY) in line, (
        f"the message does not name {MEMORY}, so a user cannot tell what to "
        f"change {INVALID} to: {line!r}")
    assert "even" not in line.lower(), (
        f"the message still says 'even', which is the criterion that was "
        f"wrong: {line!r}")


def test_a_multiple_of_the_block_size_still_runs(accepted):
    assert accepted.returncode == 0, (
        f"{VALID} moments is a multiple of {MEMORY} and must be accepted; the "
        f"check must not have become stricter than the loop requires")
    assert "Done" in accepted.stdout
