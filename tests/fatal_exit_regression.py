"""An invalid input must fail the same way every time, and say why.

KITEx validates some of its input *inside* the OpenMP region. A guard there
cannot use `exit()`: that runs the static destructors, tearing down HDF5 while
the sibling threads are still in it. The process then dies with SIGSEGV instead
of the code the guard asked for -- and the diagnostic, still sitting in `cout`'s
buffer, is never printed at all, so the user gets a bare segfault with no clue
what was wrong. Measured before the fix: five runs of one invalid file gave
139, 139, 139, 139, 139 and printed nothing past "Starting program...".

`Src/Generic.hpp`'s `fatal_exit()` flushes and then calls `_Exit`, which skips
the teardown. Guards in `main.cpp` run before the region and still use `exit()`.

Reaching these guards needs an HDF5 written by hand, because the Python exporter
rejects the same input earlier with a better message -- so the test patches a
valid file, which is what any hand-authored input amounts to.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

RUNS = 5

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
lat.add_sublattices(['A', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'A', 'A', -1.0), ([0, 1], 'A', 'A', -1.0))

# divisions=[2,2] so there are sibling threads for the teardown to race with.
cfg = kite.Configuration(divisions=[2, 2], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=False,
                         precision=1, spectrum_range=[-4.5, 4.5])
calc = kite.Calculation(cfg)
calc.conductivity_dc(direction='xy', num_points=100, num_moments=192,
                     num_random=1, num_disorder=1, temperature=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def invalid_input(tmp_path_factory):
    h5py = pytest.importorskip("h5py")
    d = tmp_path_factory.mktemp("fatal")
    path = d / "c.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, str(path)], cwd=d,
                   check=True, capture_output=True, text=True)
    # Odd NumMoments trips Gamma2D's guard, which runs inside the region.
    with h5py.File(path, "r+") as f:
        f["/Calculation/conductivity_dc/NumMoments"][...] = 201
    return path


def test_the_same_invalid_input_always_fails_the_same_way(invalid_input):
    codes = []
    for _ in range(RUNS):
        codes.append(subprocess.run([str(KITEX), str(invalid_input)],
                                    capture_output=True).returncode)
    assert len(set(codes)) == 1, (
        f"the same invalid input gave exit codes {codes}. A guard inside the "
        "OpenMP region is racing the HDF5 teardown")


def test_it_fails_with_a_status_not_a_signal(invalid_input):
    """subprocess reports a signal as a negative returncode."""
    code = subprocess.run([str(KITEX), str(invalid_input)],
                          capture_output=True).returncode
    assert code > 0, (
        f"KITEx died from signal {-code} instead of exiting. exit() from inside "
        "the OpenMP region tears down HDF5 under the sibling threads")
    assert code == 1, f"expected exit 1 from the guard, got {code}"


def test_the_diagnostic_actually_reaches_the_user(invalid_input):
    """The whole point: an abort that explains itself."""
    proc = subprocess.run([str(KITEX), str(invalid_input)], capture_output=True,
                          text=True)
    assert "moments" in (proc.stdout + proc.stderr).lower(), (
        "the guard's message never appeared. It was still buffered in cout when "
        "the process went down, which is what made the old failure so opaque: "
        f"stdout was {proc.stdout.strip()!r}")
