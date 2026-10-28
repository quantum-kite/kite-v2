"""KITE-tools must decode the ldos_map site index, and decode it correctly.

`/Calculation/ldos_map/Map` is a flat array over sites; turning a column back
into (x, y[, z], orbital) is the whole job of this reader, and getting it wrong
produces a plausible map with the values on the wrong atoms -- the worst kind of
failure for a real-space figure, because nothing looks broken.

The reference is a structural fact rather than a numerical one: a vacancy
removes the orbital, so the LDoS there is *exactly* zero. Placing one at a known
(cell, sublattice) on a lattice whose two sides differ pins the convention
uniquely -- `x + y*Lx + orbital*Nt`, x fastest and orbital slowest. Every other
candidate ordering maps that vacancy to a different column.

The operator variant is checked against the plain map it should agree with: the
projector on sublattice A, summed over orbitals, is the A-orbital LDoS. They are
separate stochastic runs, and the random-vector stream is not seeded (`seed_h`
does not reach it), so they agree only to the sampling error and every run draws
differently. The comparison is therefore a z-score against the standard errors
KITEx stores alongside the values, which is scale-free and does not need a
tolerance picked by hand -- an earlier version compared the means with a 5%
tolerance and flaked, because `l1`'s mean is two orders of magnitude smaller than
`l0`'s. As a bonus this also checks the quoted errors are honest: if they were
understated the z-scores would not sit near 1.
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

LX, LY = 16, 24          # multiples of TILE=8, and different, so no ordering ties
VACANCY_CELL = (5, 3)    # also different from each other
NT = LX * LY

_BUILDER = r'''
import sys, numpy as np, kite
from kite import lattice as latt

LX, LY = 16, 24
lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
lat.add_sublattices(['A', [0, 0], +0.3], ['B', [0.5, 0.5], -0.3])
lat.add_hoppings(([0, 0], 'A', 'B', -1.0), ([-1, 0], 'A', 'B', -1.0),
                 ([0, -1], 'A', 'B', -1.0), ([-1, -1], 'A', 'B', -1.0))

dis = kite.StructuralDisorder(lat, position=[[5, 3]])
dis.add_vacancy('A')

cfg = kite.Configuration(divisions=[1, 1], length=[LX, LY],
                         boundaries=['periodic', 'periodic'], is_complex=True,
                         precision=1, spectrum_range=[-5, 5])
calc = kite.Calculation(cfg)
calc.add_orbital_index('A', 0)
calc.add_orbital_index('B', 1)
calc.add_orbital_coupling('A', 'A', 1.0, 'l0')
calc.add_orbital_coupling('B', 'B', 1.0, 'l1')
calc.ldos_map(energy_=0.2, sigma_=0.15, vectors_=96, operators=['l0', 'l1'])
kite.config_system(lat, cfg, calc, filename=sys.argv[1],
                   disorder_structural=[dis])
'''


def _columns(path):
    rows = [line.split() for line in path.read_text().splitlines()
            if not line.startswith("#")]
    assert rows, f"{path.name} has no data rows"
    return rows


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("ldosmap")
    subprocess.run([sys.executable, "-c", _BUILDER, "m.h5"], cwd=d,
                   check=True, capture_output=True, text=True)
    subprocess.run([str(KITEX), "m.h5"], cwd=d, check=True, capture_output=True)
    subprocess.run([str(KITE_TOOLS), "m.h5", "--LDOSmap", "-X", "-N", "m.dat"],
                   cwd=d, check=True, capture_output=True)
    return d


def test_the_vacancy_site_is_decoded_to_the_right_coordinates(run):
    rows = _columns(run / "m.dat")
    assert len(rows) == NT * 2, f"expected one row per orbital per cell, got {len(rows)}"

    zeros = [(int(r[0]), int(r[1]), int(r[2])) for r in rows if float(r[3]) == 0.0]
    assert zeros == [(VACANCY_CELL[0], VACANCY_CELL[1], 0)], (
        f"the exactly-zero site decoded to {zeros}, not the vacancy at "
        f"{VACANCY_CELL} on orbital 0. The flat index is x + y*Lx + orbital*Nt; "
        "another ordering puts the map on the wrong atoms")


def test_the_other_orbital_of_the_vacancy_cell_is_not_zero(run):
    """Guards against a decode that zeroes a whole cell rather than one orbital."""
    rows = _columns(run / "m.dat")
    partner = [float(r[3]) for r in rows
               if (int(r[0]), int(r[1]), int(r[2])) == (*VACANCY_CELL, 1)]
    assert partner and partner[0] > 0.0, (
        "the B orbital of the vacancy cell is also zero; the vacancy removes one "
        "orbital, not the cell")


def test_the_energy_is_written_back_in_eV(run):
    """Energy and width are stored rescaled; the header must undo that."""
    header = (run / "m.dat").read_text().splitlines()[0]
    match = re.search(r"at E = ([-\d.eE+]+) eV, sigma = ([-\d.eE+]+) eV", header)
    assert match, f"the header no longer states E and sigma in eV: {header!r}"
    assert float(match.group(1)) == pytest.approx(0.2, abs=1e-9)
    assert float(match.group(2)) == pytest.approx(0.15, abs=1e-9)


def test_the_units_are_stated(run):
    text = (run / "m.dat").read_text()
    assert "1/eV" in text, "the gaussian branch no longer states its units"


def test_the_operator_map_agrees_with_the_plain_map(run):
    """A sublattice projector reproduces that sublattice's LDoS."""
    rows = _columns(run / "m.dat")
    orbital = np.array([int(r[2]) for r in rows])
    value = np.array([float(r[3]) for r in rows])
    error = np.array([float(r[4]) for r in rows])

    for label, which in (("l0", 0), ("l1", 1)):
        path = run / f"m_{label}.dat"
        assert path.exists(), f"--LDOSmap wrote no file for operator {label}"
        op_rows = _columns(path)
        assert len(op_rows) == NT, (
            f"{path.name} has {len(op_rows)} rows; the operator map carries one "
            "value per unit cell, summed over orbitals")
        op = np.array([float(r[3]) for r in op_rows])
        op_err = np.array([float(r[4]) for r in op_rows])
        plain = value[orbital == which]
        plain_err = error[orbital == which]

        # Site by site, in units of the combined standard error. For two
        # unbiased estimates of the same quantity this is a standard normal, so
        # its RMS is 1 and its median absolute value is 0.67.
        # The vacancy site carries value and error exactly zero in both maps,
        # so it has no z-score. Drop it, and check nothing else was dropped.
        spread = np.sqrt(op_err ** 2 + plain_err ** 2)
        live = spread > 0
        assert (~live).sum() == (1 if which == 0 else 0), (
            f"{(~live).sum()} sites have zero error in the {label} comparison; "
            "only the vacancy on orbital 0 should")
        z = (op[live] - plain[live]) / spread[live]
        assert np.sqrt((z ** 2).mean()) < 3.0, (
            f"the {label} projector map disagrees with the orbital-{which} LDoS "
            f"by {np.sqrt((z ** 2).mean()):.1f} standard errors RMS; they "
            "estimate the same quantity")
        assert np.median(np.abs(z)) < 2.0, (
            f"the {label} projector map is systematically offset from the "
            f"orbital-{which} LDoS (median |z| = {np.median(np.abs(z)):.1f})")

    # And the two sublattices must not be interchangeable, or the comparison
    # above would pass with the operators swapped.
    a = value[orbital == 0].mean()
    b = value[orbital == 1].mean()
    assert abs(a - b) / max(a, b) > 0.5, (
        "the two sublattices carry nearly the same weight, so this test cannot "
        "tell a correct operator mapping from a swapped one")
