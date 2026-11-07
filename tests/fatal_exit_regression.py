"""An invalid input must fail the same way every time, and say why.

KITEx validates some of its input *inside* the OpenMP region. A guard there
cannot use `exit()`: that runs the static destructors, tearing down HDF5 while
the sibling threads are still in it. The process then dies with SIGSEGV instead
of the code the guard asked for -- and the diagnostic, still sitting in `cout`'s
buffer, is never printed at all, so the user gets a bare segfault with no clue
what was wrong. Measured before the fix: five runs of one invalid file gave
139, 139, 139, 139, 139 and printed nothing past "Starting program...".

A `fatal_exit()` that flushed and then called `_Exit`, skipping the teardown,
existed in `Src/Generic.hpp` and was removed when the sources were restored. The
guards in `Gamma2D.cpp` use `exit()` again, so the race is live: measured
2026-11-06, twelve runs of one invalid file died of SIGSEGV eleven times with
nothing printed, and exited 1 once with the message. Run from a shell that shows
as 139; `subprocess` reports the same thing as -11.

The test patches a valid file rather than writing one by hand, but that is
convenience, not necessity: `num_moments=201` passes straight through
`Calculation.conductivity_dc` and the exporter, so ordinary Python reaches this.

Each case therefore runs KITEx `RUNS` times and requires the right answer every
time. Asking only that the runs *agree* is not enough -- at roughly eleven in
twelve they agree on the wrong answer, and the case passes while reporting a
determinism that is not there. That is what made this file look like a flaky
test rather than a detector.
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


def _codes(path, runs=RUNS):
    return [subprocess.run([str(KITEX), str(path)],
                           capture_output=True).returncode for _ in range(runs)]


def test_the_same_invalid_input_always_fails_the_same_way(invalid_input):
    """Agreement is necessary but not sufficient -- they must agree on 1.

    Checking only `len(set(codes)) == 1` passes whenever every run happens to
    land on the same side of the race, which is the common case.
    """
    codes = _codes(invalid_input)
    assert set(codes) == {1}, (
        f"the same invalid input gave exit codes {codes}; expected 1 every "
        "time. A negative code is a fatal signal -- -11 is SIGSEGV, which a "
        "shell reports as 139. A guard inside the OpenMP region called exit(), "
        "and its static destructors tore down HDF5 under the sibling threads. "
        "A mixture of 1 and -11 is the race caught in the act; all -11 is it "
        "losing every time")


def test_it_fails_with_a_status_not_a_signal(invalid_input):
    """subprocess returns -signum for a fatal signal; a shell shows 128+signum."""
    codes = _codes(invalid_input)
    crashed = [c for c in codes if c < 0 or c == 139]
    assert not crashed, (
        f"{len(crashed)} of {len(codes)} runs died from a signal rather than "
        f"exiting: {codes}. The guard is reached inside the parallel region and "
        "exit() runs the static destructors there")
    assert set(codes) == {1}, f"expected exit 1 from the guard, got {codes}"


def test_the_diagnostic_actually_reaches_the_user(invalid_input):
    """The whole point: an abort that explains itself, every time."""
    missing = []
    for _ in range(RUNS):
        proc = subprocess.run([str(KITEX), str(invalid_input)],
                              capture_output=True, text=True)
        if "moments" not in (proc.stdout + proc.stderr).lower():
            missing.append(proc.returncode)
    assert not missing, (
        f"{len(missing)} of {RUNS} runs printed no diagnostic, exiting "
        f"{missing}. The guard's message was still buffered in cout when the "
        "process went down -- which is what makes this failure opaque to a user: "
        "an invalid input produces a bare crash with no reason given")
