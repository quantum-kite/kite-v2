"""`Gamma1D` and ARPES must refuse an odd moment count.

Both advance two moments per iteration and write both:

    for (int m = 0; m < N_moments; m += 2)
        gamma.block(0, m, 1, 2) += ...      // writes m and m+1

so an odd count makes the last pass write one element past the end. Neither
checked. This is the converse of the `Gamma2D` defect: there the guard existed
and tested the wrong thing, here there was no guard at all.

Demonstrated before the fix by instrumenting the write: at 63 moments `Gamma1D`
writes columns 62 and 63 of an array with 63, and ARPES writes row 63 of an
array with 63 -- both with the run exiting 0. One element, not the 1607 of the
block-stepping case, which is why it never surfaced as a crash.

`Gamma1D` carries `dos`, `conductivity_optical` and the first-order part of
`conductivity_optical_nonlinear`. Python refuses an odd count for all of these,
so reaching the guard means a hand-written HDF5 -- which is what it is for.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

EVEN = 64
ODD = 63

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

which = sys.argv[2]
lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0], ['B', [0.5, 0.5], 0.0])
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([-1, 0], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0), ([-1, -1], 'A', 'B', -1.0))
cfg = kite.Configuration(divisions=[1, 1], length=[16, 16],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-5, 5],
                         seed_h=7, seed_v=7)
calc = kite.Calculation(cfg)
if which == 'dos':
    calc.dos(num_points=256, num_moments=64, num_random=1, num_disorder=1)
else:
    calc.arpes(k_vector=np.array([[0.0, 0.0], [0.2, 0.2]]), weight=[1.0, 1.0],
               num_moments=64, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''

#: which calculation -> the HDF5 path holding its moment count
GROUP = {"dos": "/Calculation/dos/NumMoments",
         "arpes": "/Calculation/arpes/NumMoments"}


def _run(tmp, which, moments):
    h5py = pytest.importorskip("h5py")
    path = tmp / f"{which}_{moments}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path), which], cwd=tmp,
                   check=True, capture_output=True, text=True)
    if moments != EVEN:
        with h5py.File(path, "r+") as f:
            f[GROUP[which]][...] = moments
    return path, subprocess.run([str(KITEX), str(path)], capture_output=True,
                                text=True)


@pytest.mark.parametrize("which", ["dos", "arpes"])
def test_an_odd_count_is_refused(tmp_path, which):
    path, proc = _run(tmp_path, which, ODD)
    assert proc.returncode == 1, (
        f"{which} accepted {ODD} moments; the last pass writes one element "
        f"past the array (exit {proc.returncode})")
    out = proc.stdout + proc.stderr
    assert "even" in out.lower() and str(ODD) in out, (
        f"the message must say what is wrong and name the value: {out.strip()!r}")


@pytest.mark.parametrize("which", ["dos", "arpes"])
def test_the_refusal_happens_before_anything_is_written(tmp_path, which):
    h5py = pytest.importorskip("h5py")
    path, _ = _run(tmp_path, which, ODD)
    dataset = {"dos": "/Calculation/dos/MU", "arpes": "/Calculation/arpes/kMU"}
    with h5py.File(path, "r") as f:
        assert dataset[which] not in f, (
            "aborting after writing the results is no better than not aborting")


@pytest.mark.parametrize("which", ["dos", "arpes"])
def test_an_even_count_still_runs(tmp_path, which):
    """A guard that over-rejects would be its own defect."""
    _, proc = _run(tmp_path, which, EVEN)
    assert proc.returncode == 0, f"{which} refused a valid count: {proc.stdout}"
