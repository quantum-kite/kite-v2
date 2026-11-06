"""A KITE-tools flag is given in energy units, and the reader rescales it.

Every reader keeps its parameters in KPM reduced units and converts on the way
in: `variables.X_Temp/scale`, `(variables.X_Fermi - shift)/scale`. A reader that
forgets the conversion still accepts the flag, still runs, and still prints a
plausible number -- it is simply answering a different question. There is no
banner that shows the discrepancy, so only a comparison finds it.

Two such defects are pinned here, both measured before they were fixed:

`--ARPES -T` assigned `variables.ARPES_Temp` with no `/scale`, so it was off by
`energy_scale` while every other reader took eV. ARPES prints no
temperature line at all, so nothing on screen revealed it.

`custom_two`'s default Fermi range was `+-1.0/scale`, omitting the `- shift`
that `conductivity_dc` applies and that `custom_two`'s own `-F` override already
applied, so on a shifted spectrum the default swept a window displaced by the
shift.

Both are checked against a reference rather than a stored number: a flag set to
a reader's own default must reproduce the default run, and two readers given the
same window must sweep the same window.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"
KITE_TOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITE_TOOLS.exists()),
    reason="needs a built KITEx and KITE-tools")


def _numeric(path):
    """These writers interleave headers with the data; keep the widest rows."""
    rows = []
    for line in Path(path).read_text().splitlines():
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            continue
    width = max(len(r) for r in rows)
    return np.array([r for r in rows if len(r) == width])


# --- --ARPES -T is in energy units -------------------------------------------

_ARPES = r'''
import sys, numpy as np, kite
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
k = np.array([[0.0, 0.0], [0.3, 0.0], [0.5, 0.5]])
calc.arpes(k_vector=k, weight=[0.5, 1.5], num_moments=32, num_disorder=1)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''

#: arpes.cpp's built-in default is `temperature = 0.1/scale`, i.e. 0.1 eV.
ARPES_DEFAULT_EV = 0.1
ARPES_SCALE = 4.1


@pytest.fixture(scope="module")
def arpes(tmp_path_factory):
    d = tmp_path_factory.mktemp("arpes_units")
    subprocess.run([sys.executable, "-c", _ARPES, "a.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "a.h5"], cwd=d, check=True, capture_output=True)

    def run(name, *extra):
        subprocess.run([str(KITE_TOOLS), "a.h5", "--ARPES", "-N", name, *extra],
                       cwd=d, check=True, capture_output=True, text=True)
        return _numeric(d / (name + ".dat"))   # this writer appends .dat
    return run


def test_the_arpes_temperature_flag_is_in_energy_units(arpes):
    """`-T` equal to the built-in default must reproduce the default run."""
    default = arpes("d")
    asked = arpes("a", "-T", str(ARPES_DEFAULT_EV))
    assert np.array_equal(default, asked), (
        f"--ARPES -T {ARPES_DEFAULT_EV} does not reproduce the built-in default "
        f"of {ARPES_DEFAULT_EV} eV, so -T is not being read in energy units")


def test_the_arpes_temperature_flag_is_not_in_reduced_units(arpes):
    """The old behaviour, which must not come back."""
    default = arpes("d2")
    reduced = arpes("r", "-T", str(ARPES_DEFAULT_EV / ARPES_SCALE))
    assert not np.array_equal(default, reduced), (
        "-T given in reduced units reproduces the default, so the reader is "
        "skipping the /energy_scale conversion again")


# --- The default Fermi window is the same for both readers -------------------

_SHIFTED = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

a, acc = 0.24595, 0.142
lat = latt.Lattice(a1=a*np.array([1.0, 0.0]), a2=a*np.array([0.5, 0.5*np.sqrt(3)]))
# A uniform onsite term puts the spectrum off centre, so energy_shift != 0.
lat.add_sublattices(('A', [0, -acc/2], 3.0), ('B', [0, acc/2], 3.0))
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([1, -1], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0))

cfg = kite.Configuration(divisions=[2, 2], length=[32, 32],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-1.0, 7.0],
                         seed_h=20261031, seed_v=20261031)
calc = kite.Calculation(cfg)
calc.conductivity_dc(num_points=100, num_moments=64, num_random=1,
                     num_disorder=1, direction='xx', temperature=0.05)
A = custom.Vertex(64, [[1.0, "vx"]])
B = custom.Vertex(64, [[1.0, "vx"]])
calc.custom_two(stream_=[A, B], num_random_=1, num_disorder_=1,
                num_points_=100, temperature_=0.05)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''

SHIFT = 3.0


@pytest.fixture(scope="module")
def shifted(tmp_path_factory):
    d = tmp_path_factory.mktemp("fermi_shift")
    subprocess.run([sys.executable, "-c", _SHIFTED, "s.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "s.h5"], cwd=d, check=True, capture_output=True)
    # One invocation computes both; custom_two writes its own default filename.
    subprocess.run([str(KITE_TOOLS), "s.h5", "--CondDC", "-N", "dc.dat"],
                   cwd=d, check=True, capture_output=True, text=True)
    return d


def test_the_shift_is_actually_present(shifted):
    """Without a shift this comparison cannot distinguish anything."""
    import h5py
    with h5py.File(shifted / "s.h5", "r") as f:
        assert float(np.array(f["EnergyShift"])) == pytest.approx(SHIFT)


def test_both_readers_sweep_the_same_default_fermi_window(shifted):
    """custom_two dropped `- shift`, so its window slid by the shift."""
    dc = _numeric(shifted / "dc.dat")
    custom = _numeric(shifted / "custom_two.dat")
    assert custom[0, 0] == pytest.approx(dc[0, 0], abs=1e-6), (
        f"custom_two starts its Fermi sweep at {custom[0, 0]:+.4f} eV where "
        f"conductivity_dc starts at {dc[0, 0]:+.4f} eV, a displacement of "
        f"{custom[0, 0] - dc[0, 0]:+.4f} eV against an energy shift of {SHIFT}")
    assert custom[-1, 0] == pytest.approx(dc[-1, 0], abs=1e-6)
