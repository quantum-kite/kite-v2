"""KITE-tools must write the wave-packet series on a physical time axis, without `Id`.

Two things this reader has to get right, and neither is arithmetic.

The time axis is stored *rescaled*: `timestep` is the argument of the Bessel
expansion of `exp(-i Htilde ttilde)`, so a step of `dt` in physical time is
`dt * EnergyScale` there. Leaving it alone would hand the user an axis that is
wrong by a factor of `EnergyScale` -- 5 in this test -- while still looking
entirely plausible. The reference is a dispersion-free closed form: a spin-1/2
in a uniform transverse field precesses as `<Sz(t)> = (1/2)cos(2Bt)` and
`<Sy(t)> = -(1/2)sin(2Bt)` exactly, for any packet width, because the Bloch
Hamiltonian factorizes as `eps(k) I + B sigma_x` at every k. That closed form
holds on the physical axis and on no other, so it pins the conversion.

`Id` must not appear. It is the amplitude at a probing point whose position
depends on how the lattice was split across threads, so the same physics gives
a different `Id` at a different thread count. AC-WAVE-06 rules that it is
hidden, not explained, and the obvious implementation of this reader -- loop
over the datasets in the group -- puts it straight back in front of users.
"""

import re
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

B = 0.1                  # transverse field, so the precession is at 2B
ENERGY_SCALE = 5.0       # spectrum_range=[-5, 5]
TIMESTEP = 0.5           # rescaled, i.e. 0.1 hbar/eV
NUM_POINTS = 200
CENTRE = 32              # mean_value, on a 64x64 lattice
WIDTH = 4.0

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

B, L = 0.1, 64
lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
lat.add_sublattices(['up', [0, 0], 0.0], ['down', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'up', 'up', -1.0), ([0, 1], 'up', 'up', -1.0),
                 ([1, 0], 'down', 'down', -1.0), ([0, 1], 'down', 'down', -1.0),
                 ([0, 0], 'up', 'down', B))

cfg = kite.Configuration(divisions=[1, 1], length=[L, L],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-5, 5])
calc = kite.Calculation(cfg)
calc.add_orbital_index('up', 0)
calc.add_orbital_index('down', 1)
calc.add_orbital_coupling('down', 'up', 0.5, 'l0')     # Sx
calc.add_orbital_coupling('up', 'down', 0.5, 'l0')
calc.add_orbital_coupling('down', 'up', -0.5j, 'l1')   # Sy
calc.add_orbital_coupling('up', 'down', 0.5j, 'l1')
calc.add_orbital_coupling('up', 'up', 0.5, 'l2')       # Sz
calc.add_orbital_coupling('down', 'down', -0.5, 'l2')

calc.gaussian_wave_packet(
    num_points=200, num_moments=128, num_disorder=1,
    k_vector=np.array([[0.0, 0.0]]), spinor=np.array([[1.0 + 0j, 0.0 + 0j]]),
    width=4.0, timestep=0.5, mean_value=[L // 2, L // 2],
    operators=['l0', 'l1', 'l2'])
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("wavepacket")
    subprocess.run([sys.executable, "-c", _BUILDER, "w.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "w.h5"], cwd=d, check=True, capture_output=True)
    subprocess.run([str(KITE_TOOLS), "w.h5", "--WavePacket", "-X", "-N", "w.dat"],
                   cwd=d, check=True, capture_output=True)
    return d / "w.dat"


@pytest.fixture(scope="module")
def data(run):
    return np.loadtxt(run)


def test_the_probing_point_amplitude_is_not_emitted(run, data):
    """`Id` depends on the thread decomposition, so it is hidden (AC-WAVE-06)."""
    header = "\n".join(line for line in run.read_text().splitlines()
                       if line.startswith("#"))
    assert not re.search(r"\bId\b", header), (
        "the column list names Id. Its probing point moves with the thread "
        "decomposition, so the same physics gives a different Id at a different "
        "thread count; it is hidden, not explained")
    # t, t_rescaled, <x>, <y>, Varx, Vary, and the three spin components.
    assert data.shape == (NUM_POINTS, 9), (
        f"expected 9 columns, got {data.shape[1]} -- an extra one is most "
        "likely Id leaking back in")


def test_the_time_axis_is_physical_not_rescaled(data):
    """The closed-form precession holds on t, and on no other axis."""
    t, t_rescaled = data[:, 0], data[:, 1]
    sz, sy = data[:, 8], data[:, 7]

    assert np.abs(sz - 0.5 * np.cos(2 * B * t)).max() < 1e-4, (
        "<Sz> does not follow (1/2)cos(2Bt) on the physical axis")
    assert np.abs(sy + 0.5 * np.sin(2 * B * t)).max() < 1e-4, (
        "<Sy> does not follow -(1/2)sin(2Bt) on the physical axis")

    # And it must not hold on the stored axis, or the conversion is untested.
    assert np.abs(sz - 0.5 * np.cos(2 * B * t_rescaled)).max() > 0.1, (
        "the two axes agree, so this test cannot tell a converted time axis "
        "from an unconverted one")
    assert t_rescaled == pytest.approx(t * ENERGY_SCALE, abs=1e-9)
    assert t_rescaled[1] == pytest.approx(TIMESTEP)


def test_the_field_free_component_stays_zero(data):
    """Sx commutes with the field, and the packet starts with none of it."""
    assert np.abs(data[:, 6]).max() < 1e-12


def test_the_packet_starts_where_it_was_seeded(data):
    """Pins the position and variance columns, and their order."""
    assert data[0, 2] == pytest.approx(CENTRE)
    assert data[0, 3] == pytest.approx(CENTRE)
    assert data[0, 4] == pytest.approx(WIDTH ** 2 / 2)
    assert data[0, 5] == pytest.approx(WIDTH ** 2 / 2)
    # k = 0, so the packet spreads without drifting. It does not return to the
    # centre to the last digit -- it has spread far enough to wrap the periodic
    # cell by the end -- so this is a drift check, not an identity.
    assert data[-1, 2] == pytest.approx(CENTRE, abs=1e-2)
    assert data[-1, 4] > 3 * data[0, 4], (
        "the variance column does not grow; the packet disperses")


def test_the_units_are_stated(run):
    text = run.read_text()
    assert "hbar/eV" in text, "the header no longer states the time unit"
    assert "t_rescaled" in text, (
        "the header no longer says that the second column is the stored axis")
    assert f"EnergyScale = {ENERGY_SCALE:g} eV" in text, (
        "the header no longer quotes the scale that relates the two axes")
