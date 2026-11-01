"""`-T` means the same thing in every reader that accepts it.

Four defects of this cycle were one defect wearing four hats: the energy
scale/shift convention is applied per reader, by hand, with nothing enforcing
that the readers agree. `--ARPES -T` was read in reduced units (BUG-013);
`custom_two`'s default Fermi window dropped the shift (BUG-014); `custom_two`
divided a temperature that Python had already divided (BUG-005); a solver
dropped the shift outright (BUG-010). Each was closed by making one file match
its siblings, and each was found by looking at that file.

So this sweeps instead of naming. Every reader that accepts `-T` is run three
ways on one exported file:

  * with no flag, where the temperature comes from the HDF5 file;
  * with `-T` set to that same temperature in eV, which must reproduce the
    no-flag run **byte for byte**;
  * with `-T` set to that temperature already divided by `energy_scale`, which
    must not -- a reader that skips its own conversion would accept the reduced
    number as if it were eV.

`--ARPES` has no `Temperature` dataset to fetch, so its no-flag run uses the
literal in `arpes.cpp`; that literal is read out of the source here rather than
written down, so moving it breaks this file instead of silently making the
comparison vacuous.

`OMP_NUM_THREADS=1` is not decoration. KITE-tools accumulates its integrals
with no fixed reduction order, so two runs of the *same* binary on the same
input differ in the last digits at more than one thread, and a byte comparison
would be measuring the scheduler.

Writing it as a sweep paid for itself immediately: `--CondOpt` failed the first
rule. It read the stored temperature into `temperature` but never recomputed
`beta`, which is the only thing its integrand uses, so every run computed at
the hardcoded 0.001 eV while printing the temperature that had been asked for.
That is BUG-015, now fixed; this file is what found it and what would find the
next one.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")

#: The temperature the exported calculations ask for, in eV. Far from the
#: 0.001 eV that every reader falls back to, so a reader using the fallback is
#: unmistakable rather than a rounding argument.
TEMPERATURE_EV = 0.5

_EXPORT = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

a, acc = 0.24595, 0.142
lat = latt.Lattice(a1=a*np.array([1.0, 0.0]), a2=a*np.array([0.5, 0.5*np.sqrt(3)]))
lat.add_sublattices(('A', [0, -acc/2], 0.0), ('B', [0, acc/2], 0.0))
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([1, -1], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0))

cfg = kite.Configuration(divisions=[1, 1], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-4.1, 4.1],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
T = float(sys.argv[2])
calc.conductivity_dc(direction='xx', num_points=50, num_moments=32,
                     num_random=1, num_disorder=1, temperature=T)
calc.conductivity_optical(direction='xx', num_points=50, num_moments=32,
                          num_random=1, num_disorder=1, temperature=T)
calc.conductivity_optical_nonlinear(direction='xxx', num_points=50, num_moments=32,
                                    num_random=1, num_disorder=1, temperature=T)
A = custom.Vertex(32, [[1.0, "vx"]])
B = custom.Vertex(32, [[1.0, "vx"]])
calc.custom_two(stream_=[A, B], num_random_=1, num_disorder_=1,
                num_points_=50, temperature_=T)
calc.arpes(k_vector=np.array([[0.0, 0.0], [0.3, 0.0], [0.5, 0.5]]),
           weight=[0.5, 1.5], num_moments=32, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _arpes_default_ev():
    """arpes.cpp has no Temperature dataset; read its literal from the source."""
    text = (REPO / "tools" / "Src" / "Spectral" / "arpes.cpp").read_text()
    match = re.search(r"temperature\s*=\s*([0-9.]+)\s*/\s*scale\s*;", text)
    assert match, "arpes.cpp no longer sets a literal default temperature"
    return float(match.group(1))


#: (flag, output file, the temperature its no-flag run actually uses, in eV)
READERS = [
    ("--CondDC", "condDC.dat", lambda: TEMPERATURE_EV),
    ("--CondOpt", "optcond.dat", lambda: TEMPERATURE_EV),
    ("--CondOpt2", "nonlinear_cond.dat", lambda: TEMPERATURE_EV),
    ("--CustomTwo", "custom_two.dat", lambda: TEMPERATURE_EV),
    ("--ARPES", "arpes.dat", _arpes_default_ev),
]


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    """One export and one KITEx run; the readers only re-read it."""
    import h5py
    directory = tmp_path_factory.mktemp("temperature_sweep")
    subprocess.run([sys.executable, "-c", _EXPORT, "t.h5", str(TEMPERATURE_EV)],
                   cwd=directory, check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "t.h5"], cwd=directory, check=True, capture_output=True)
    with h5py.File(directory / "t.h5", "r") as handle:
        scale = float(handle["EnergyScale"][()])
    return directory, scale


def _read(exported, flag, tag, *extra):
    directory, _ = exported
    work = directory / f"{flag.lstrip('-')}_{tag}"
    work.mkdir()
    shutil.copy(directory / "t.h5", work / "t.h5")
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    result = subprocess.run([str(KITE_TOOLS), "t.h5", flag, *extra], cwd=work,
                            env=environment, capture_output=True, text=True)
    assert result.returncode == 0, f"{flag} {extra} failed: {result.stdout[-500:]}"
    return work


def test_the_export_carries_every_reader_this_file_sweeps(exported):
    """Otherwise a reader could drop out and its cases would pass vacuously."""
    import h5py
    directory, _ = exported
    with h5py.File(directory / "t.h5", "r") as handle:
        present = set(handle["/Calculation"].keys())
    assert {"conductivity_dc", "conductivity_optical",
            "conductivity_optical_nonlinear", "CustomTwo", "arpes"} <= present, (
        f"the export is missing calculations this file sweeps: {sorted(present)}")


@pytest.mark.parametrize("flag,output,temperature", READERS)
def test_the_flag_at_its_own_value_reproduces_the_default_run(flag, output,
                                                              temperature, exported):
    asked = temperature()
    plain = _read(exported, flag, "plain")
    flagged = _read(exported, flag, "same", "-T", repr(asked))
    assert (plain / output).read_bytes() == (flagged / output).read_bytes(), (
        f"{flag} -T {asked} does not reproduce the run with no flag, although "
        f"that run uses {asked} eV. One of the two paths is not applying the "
        "energy-scale conversion the other applies")


@pytest.mark.parametrize("flag,output,temperature", READERS)
def test_the_flag_is_not_read_in_reduced_units(flag, output, temperature, exported):
    """The shape BUG-013 had: the flag taken as already divided."""
    _, scale = exported
    reduced = temperature() / scale
    plain = _read(exported, flag, "plain2")
    flagged = _read(exported, flag, "reduced", "-T", repr(reduced))
    assert (plain / output).read_bytes() != (flagged / output).read_bytes(), (
        f"{flag} -T {reduced} reproduces the default run, so the flag is being "
        "read in reduced units instead of eV")
