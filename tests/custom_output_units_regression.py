"""Every custom_* result must say what its number means, or say it is uncalibrated.

`custom_one.dat` and `custom_two.dat` used to be bare columns of numbers. That is
the worst case for a quantity whose absolute scale depends on what the user put in
the vertex: the file looks finished, and nothing in it says that a vertex built
from KITE's raw `v` token is a commutator with the *rescaled* Hamiltonian and so
needs a factor of `EnergyScale` per token before it means anything physical. One
shipped processor already supplies that factor by hand
(`examples/process_haldane_orbital_magnetization.py`); nothing told the next user to.

The `custom_two` case is the opposite problem -- it IS calibrated, at
`NumVelocities = 2`, but only with a convention the file never named. A spin
conductivity quoted in $e^2/h$ has $\\hbar/2e$ absorbed, so the Kane-Mele plateau
reads 2 and not 1, and a reader comparing against a paper has no way to know which
of the two they are holding.

These assertions are about what the header says, not about the numbers; the numbers
have their own tests. What they protect is the pairing: a header that stops matching
the code is worse than no header.
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

#: `vertex` is spliced into the custom_one Vertex below. The first carries no
#: velocity token, the second the symmetrized spin current (1/2){v_x, s_z}.
VERTICES = {
    "plain": "[[1.0, 'l0']]",
    "current": "[[0.5, 'vx.l0'], [0.5, 'l0.vx']]",
}

_BUILDER = r'''
import sys, numpy as np, kite
from kite import custom
from kite import lattice as latt

lat = latt.Lattice(a1=[1.0, 0], a2=[0, 1.0])
lat.add_sublattices(['up', [0, 0], 0.0], ['dn', [0, 0], 0.0])
lat.add_hoppings(([1, 0], 'up', 'up', -1.0), ([0, 1], 'up', 'up', -1.0),
                 ([1, 0], 'dn', 'dn', -1.0), ([0, 1], 'dn', 'dn', -1.0),
                 ([0, 0], 'up', 'dn', 0.2))

cfg = kite.Configuration(divisions=[1, 1], length=[16, 16],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-6, 6])
calc = kite.Calculation(cfg)
calc.add_orbital_index('up', 0)
calc.add_orbital_index('dn', 1)
calc.add_orbital_coupling('up', 'up', 0.5, 'l0')     # s_z = sigma_z/2
calc.add_orbital_coupling('dn', 'dn', -0.5, 'l0')

M = 64
calc.custom_one(stream_=custom.Vertex(M, VERTEX), num_random_=1, num_disorder_=1)
calc.custom_one_local(stream_=custom.Vertex(M, [[1.0, 'l0']]),
                      energy_=np.array([0.0]), position_=[[8, 8]],
                      sublattice_=['up'])
calc.custom_two(stream_=[custom.Vertex(M, [[0.5, 'vx.l0'], [0.5, 'l0.vx']]),
                         custom.Vertex(M, [[1.0, 'vy']])],
                num_random_=1, num_disorder_=1, num_points_=200, temperature_=0.01)
kite.config_system(lat, cfg, calc, filename=sys.argv[1])
'''


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    out = {}
    for name, vertex in VERTICES.items():
        d = tmp_path_factory.mktemp(f"customunits_{name}")
        subprocess.run([sys.executable, "-c",
                        _BUILDER.replace("VERTEX", vertex), "c.h5"],
                       cwd=d, check=True, capture_output=True, text=True)
        subprocess.run([str(KITEX), "c.h5"], cwd=d, check=True, capture_output=True)
        subprocess.run([str(KITE_TOOLS), "c.h5"], cwd=d, check=True, capture_output=True)
        out[name] = d
    return out


def _header(path):
    return "\n".join(line for line in path.read_text().splitlines()
                     if line.startswith("#"))


FILES = ["custom_one.dat", "custom_one_local.dat", "custom_two.dat"]


@pytest.mark.parametrize("name", FILES)
def test_every_custom_output_states_its_units(runs, name):
    path = runs["plain"] / name
    assert path.exists(), f"KITE-tools wrote no {name}"
    header = _header(path)
    assert "Units" in header or "units" in header, (
        f"{name} carries no units statement. Its absolute scale depends on the "
        "vertex, so a bare column of numbers cannot be compared with anything")
    assert "col 1" in header, f"{name} does not say what its columns are"


@pytest.mark.parametrize("name", FILES)
def test_the_headers_do_not_break_the_consumers(runs, name):
    """Every shipped processor reads these with np.loadtxt, which skips '#'."""
    data = np.loadtxt(runs["plain"] / name)
    assert data.ndim == 2 and len(data) > 0


def test_the_velocity_rescaling_is_stated_only_when_there_is_one(runs):
    """custom_one applies no EnergyScale power; the file has to say so."""
    with_v = _header(runs["current"] / "custom_one.dat")
    assert "EnergyScale^1" in with_v, (
        "the vertex carries one raw v token, which is a commutator with the "
        "RESCALED Hamiltonian, and custom_one applies no EnergyScale power of its "
        "own -- unlike custom_two. The file must say to multiply it back in")

    without_v = _header(runs["plain"] / "custom_one.dat")
    assert "EnergyScale" not in without_v, (
        "a vertex with no velocity token needs no rescaling, so saying to "
        "multiply by EnergyScale^0 is noise")


def test_the_spin_conductivity_convention_is_named(runs):
    """e^2/h with hbar/2e absorbed, so the quantized plateau reads 2, not 1."""
    header = _header(runs["plain"] / "custom_two.dat")
    assert "e^2/h" in header, "the custom_two header no longer states its units"
    assert "e/4pi" in header and "hbar/2e" in header, (
        "the header names e^2/h but not which convention absorbs hbar/2e, so a "
        "reader cannot tell whether a quantized plateau should be 1 or 2")
