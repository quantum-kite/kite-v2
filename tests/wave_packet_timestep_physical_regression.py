"""`timestep_physical` must mean physical time, whatever the energy scale is.

`timestep` is in rescaled units because that is what KITEx feeds to the Bessel
expansion of `exp(-i Htilde ttilde)`, and it stays that way in the file. The
cost is that the same `timestep` is a different physical step under a different
`spectrum_range`, and under an estimated one the user cannot know the step at
all until after the run. `timestep_physical` closes that: Python multiplies by
the final energy scale on export, after an unset `spectrum_range` has been
estimated.

The reference is the same dispersion-free closed form the rescaled route is
pinned to -- a spin-1/2 in a uniform transverse field precesses as
`<Sz(t)> = (1/2)cos(2Bt)` exactly, for any packet width. The test that matters
is the one with no `spectrum_range` at all: there the rescaled step is a number
nobody wrote down, and Larmor still has to come out right on the axis the user
asked for.
"""

import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

B = 0.1
PHYSICAL_STEP = 0.1
NUM_POINTS = 200

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

B, L = 0.1, 64
step_kind, step_value, scale = sys.argv[2], float(sys.argv[3]), sys.argv[4]

lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
lat.add_sublattices(['up', [0, 0], 0.0], ['down', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'up', 'up', -1.0), ([0, 1], 'up', 'up', -1.0),
                 ([1, 0], 'down', 'down', -1.0), ([0, 1], 'down', 'down', -1.0),
                 ([0, 0], 'up', 'down', B))

extra = {} if scale == 'estimated' else {'spectrum_range': [-float(scale), float(scale)]}
cfg = kite.Configuration(divisions=[1, 1], length=[L, L],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, **extra)
calc = kite.Calculation(cfg)
calc.add_orbital_index('up', 0)
calc.add_orbital_index('down', 1)
calc.add_orbital_coupling('down', 'up', 0.5, 'l0')
calc.add_orbital_coupling('up', 'down', 0.5, 'l0')
calc.add_orbital_coupling('down', 'up', -0.5j, 'l1')
calc.add_orbital_coupling('up', 'down', 0.5j, 'l1')
calc.add_orbital_coupling('up', 'up', 0.5, 'l2')
calc.add_orbital_coupling('down', 'down', -0.5, 'l2')

calc.gaussian_wave_packet(
    num_points=200, num_moments=128, num_disorder=1,
    k_vector=np.array([[0.0, 0.0]]), spinor=np.array([[1.0 + 0j, 0.0 + 0j]]),
    width=4.0, mean_value=[L // 2, L // 2], operators=['l0', 'l1', 'l2'],
    **{step_kind: step_value})
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


def _build(d, name, step_kind, step_value, scale):
    subprocess.run([sys.executable, "-c", _BUILDER, name, step_kind,
                    str(step_value), scale],
                   cwd=d, check=True, capture_output=True, text=True)
    return d / name


def _read(path):
    with h5py.File(path, "r") as f:
        g = f["Calculation"]["gaussian_wave_packet"]
        out = {"scale": float(np.array(f["EnergyScale"])),
               "timestep": float(np.array(g["timestep"]).ravel()[0]),
               "timestep_physical": float(np.array(g["timestep_physical"]).ravel()[0])}
        if "l2" in g:
            out["sz"] = np.array(g["l2"]).flatten().real
    return out


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    d = tmp_path_factory.mktemp("timestep_physical")
    out = {}
    for tag, kind, value, scale in (
            ("fixed_physical", "timestep_physical", PHYSICAL_STEP, "5"),
            ("fixed_rescaled", "timestep", PHYSICAL_STEP * 5, "5"),
            ("estimated_physical", "timestep_physical", PHYSICAL_STEP, "estimated")):
        h5 = _build(d, tag + ".h5", kind, value, scale)
        subprocess.run([str(KITEX), h5.name], cwd=d, check=True, capture_output=True)
        out[tag] = _read(h5)
    return out


def test_a_physical_step_is_rescaled_on_the_way_in(runs):
    """At a known scale the two routes have to produce the same stored step."""
    assert runs["fixed_physical"]["timestep"] == pytest.approx(
        runs["fixed_rescaled"]["timestep"]), (
        "timestep_physical=0.1 at EnergyScale=5 did not store the same rescaled "
        "step as timestep=0.5; the conversion is not timestep_physical*EnergyScale")


def test_the_rescaled_route_also_records_the_physical_step(runs):
    """One name, one convention -- both fields are written either way round."""
    assert runs["fixed_rescaled"]["timestep_physical"] == pytest.approx(PHYSICAL_STEP)


def test_an_estimated_scale_still_honours_the_physical_step(runs):
    """The case the rescaled convention cannot serve: nobody set the scale."""
    r = runs["estimated_physical"]
    assert r["scale"] != pytest.approx(5.0), (
        "the estimated bounds happen to give EnergyScale=5, so this test no "
        "longer distinguishes a converted step from an unconverted one")
    assert r["timestep_physical"] == pytest.approx(PHYSICAL_STEP)
    assert r["timestep"] == pytest.approx(PHYSICAL_STEP * r["scale"], rel=1e-6)


def test_the_precession_follows_larmor_on_the_axis_that_was_asked_for(runs):
    """AC: a user-supplied physical timestep reproduces the analytic result."""
    for tag in ("fixed_physical", "estimated_physical"):
        sz = runs[tag]["sz"]
        t = np.arange(len(sz)) * runs[tag]["timestep_physical"]
        assert np.abs(sz - 0.5 * np.cos(2 * B * t)).max() < 1e-4, (
            f"<Sz> in the {tag} run does not follow (1/2)cos(2Bt) on the "
            "physical axis the user specified")


def test_giving_both_steps_is_refused():
    """Two names for one quantity, disagreeing, must not pick a winner."""
    with pytest.raises(ValueError, match="exactly one"):
        _declare(timestep=0.5, timestep_physical=0.1)


def test_giving_neither_step_is_refused():
    """There is no default: a wave packet without a step is not a request."""
    with pytest.raises(ValueError, match="exactly one"):
        _declare()


def _declare(**step):
    import kite
    from kite import lattice as latt
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(['a', [0, 0], 0.0])
    lat.add_hoppings(([1, 0], 'a', 'a', -1.0), ([0, 1], 'a', 'a', -1.0))
    cfg = kite.Configuration(divisions=[1, 1], length=[8, 8],
                             boundaries=['periodic', 'periodic'],
                             is_complex=True, precision=1,
                             spectrum_range=[-5, 5])
    kite.Calculation(cfg).gaussian_wave_packet(
        num_points=10, num_moments=16, num_disorder=1,
        k_vector=np.array([[0.0, 0.0]]), spinor=np.array([[1.0 + 0j]]),
        width=1.0, mean_value=[4, 4], **step)
