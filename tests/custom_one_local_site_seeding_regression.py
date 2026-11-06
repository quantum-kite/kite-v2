"""`custom_one_local` must seed the site, orbital and vertex it was asked for.

It did not. `SimulationRankOneLocal.cpp` built the seed index as
`positions_[site] * r.Ld[0] + orbitals_[site]`, but `build_site` decodes its
argument in `r.Lt` coordinates, orbital-slowest, so the orbital stride is the
global cell count `r.Nt`. `r.Ld[0]` is the ghost-padded *local* extent of
direction 0 -- it appears legitimately a few lines below as the inner-product
block stride, which is likely how it got into the index.

Two distinct symptoms, and the first version of this test only caught the milder:

* at the origin, the wrong stride lands on the neighbouring cell's first
  orbital, so on a clean lattice sublattice B silently returns A's moments;
* at any other position the index runs past the array, no thread claims the
  site, and every moment comes back zero -- including mu_0, which is 1 by
  normalisation alone.

The reference is `ldos`, which accumulates `<j|T_n(H)|j>` by a verified route.
With an identity vertex `custom_one_local` computes the same object and must
agree digit for digit; with `s_z = diag(+1/2, -1/2)` it must give exactly
`+-1/2` times it, which also proves the vertex is assembled and applied
correctly rather than silently collapsing to something symmetric.

The lattice has to break the A/B symmetry or the two sublattices agree
trivially and this measures nothing: on-site +1 and -1 with
`spectrum_range=[-8, 8]` makes mu_1 and mu_3 odd in the on-site sign.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")

#: (dimension, divisions, length, position). The non-origin positions are the
#: ones that used to return all zeros; the decomposed cases put the probe in
#: different thread domains.
GEOMETRIES = [
    (2, [1, 1], [16, 16], [0, 0]),
    (2, [1, 1], [16, 16], [11, 13]),
    (2, [1, 1], [16, 16], [15, 15]),
    (2, [2, 2], [16, 16], [0, 0]),
    (2, [2, 2], [16, 16], [11, 13]),
    (2, [2, 2], [16, 16], [2, 9]),
    (3, [1, 1, 1], [16, 16, 16], [0, 0, 0]),
    (3, [2, 2, 2], [16, 16, 16], [11, 3, 13]),
]

_BUILDER = r'''
import json, sys, numpy as np, kite
from kite import lattice as latt
from kite import custom
spec = json.loads(sys.argv[1]); route = sys.argv[2]
D, div, L, pos, sub, vert = (spec[k] for k in ("D", "div", "L", "pos", "sub", "vert"))
if D == 2:
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(['A', [0, 0], +1.0], ['B', [0.5, 0.5], -1.0])
    bonds = ([0, 0], [-1, 0], [0, -1], [-1, -1])
else:
    lat = latt.Lattice(a1=[1, 0, 0], a2=[0, 1, 0], a3=[0, 0, 1])
    lat.add_sublattices(['A', [0, 0, 0], +1.0], ['B', [0.5, 0.5, 0.5], -1.0])
    bonds = ([0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, -1])
for ri in bonds:
    lat.add_hoppings([ri, 'A', 'B', -1.0])
cfg = kite.Configuration(divisions=div, length=L, boundaries=["periodic"] * D,
                         is_complex=True, precision=1, spectrum_range=[-8, 8])
calc = kite.Calculation(cfg)
energy = np.linspace(-0.1, 0.1, 3)
if route == "ldos":
    calc.ldos(energy=energy, num_moments=8, num_disorder=1,
              position=[pos], sublattice=[sub])
else:
    calc.add_orbital_index('A', 0)
    calc.add_orbital_index('B', 1)
    weights = (1.0, 1.0) if vert == "id" else (+0.5, -0.5)
    calc.add_orbital_coupling('A', 'A', weights[0], 'l0')
    calc.add_orbital_coupling('B', 'B', weights[1], 'l0')
    calc.custom_one_local(stream_=custom.Vertex(8, [[1.0, "l0"]]), energy_=energy,
                          num_disorder_=1, position_=[pos], sublattice_=[sub])
kite.config_system(lat, cfg, calc, filename=sys.argv[3])
'''


def _moments(tmp_path, spec, route):
    import h5py

    path = tmp_path / f"{route}.h5"
    subprocess.run([sys.executable, "-c", _BUILDER, json.dumps(spec), route, str(path)],
                   cwd=str(tmp_path), capture_output=True, check=True)
    subprocess.run([str(KITEX), path.name], cwd=str(tmp_path),
                   capture_output=True, check=True)
    dataset = ("/Calculation/ldos/lMU" if route == "ldos"
               else "/Calculation/CustomOneLocal/Gamma")
    with h5py.File(path, "r") as handle:
        return np.real(np.asarray(handle[dataset])).ravel()[:4]


@pytest.mark.parametrize("dimension,divisions,length,position", GEOMETRIES,
                         ids=lambda v: str(v).replace(" ", ""))
@pytest.mark.parametrize("sublattice", ["A", "B"])
@pytest.mark.parametrize("vertex", ["id", "sz"])
def test_custom_one_local_matches_ldos(tmp_path, dimension, divisions, length,
                                       position, sublattice, vertex):
    if dimension == 3 and vertex == "sz":
        pytest.skip("the spin vertex adds nothing in 3D that 2D does not cover")
    spec = dict(D=dimension, div=divisions, L=length, pos=position,
                sub=sublattice, vert=vertex)

    reference = _moments(tmp_path, spec, "ldos")
    measured = _moments(tmp_path, spec, "custom")
    expected = reference if vertex == "id" else (
        (0.5 if sublattice == "A" else -0.5) * reference)

    assert reference[0] == pytest.approx(1.0, abs=1e-9), (
        f"the ldos reference is itself wrong: mu_0 = {reference[0]}, not 1")
    assert measured == pytest.approx(expected, abs=1e-9), (
        f"custom_one_local at {position} on sublattice {sublattice} with the "
        f"{vertex} vertex gives {measured}, not {expected}")


def test_a_misseeded_probe_would_be_visible(tmp_path):
    """Guard the guard: the two sublattices must actually differ.

    If the lattice ever stops breaking the A/B symmetry, every assertion above
    passes while measuring nothing.
    """
    spec = dict(D=2, div=[1, 1], L=[16, 16], pos=[0, 0], vert="id")
    a = _moments(tmp_path, dict(spec, sub="A"), "ldos")
    b = _moments(tmp_path, dict(spec, sub="B"), "ldos")
    assert not np.allclose(a, b), (
        "the two sublattices give identical ldos moments, so this file cannot "
        "detect a probe seeded on the wrong orbital")
