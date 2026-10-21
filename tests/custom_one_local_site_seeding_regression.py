"""`custom_one_local` must seed the site and orbital it was asked for.

It does not. `Src/Simulation/Custom/SimulationRankOneLocal.cpp:147-148` builds the
seed index as `positions_[site] * r.Ld[0] + orbitals_[site]`, where `r.Ld[0]` is
the ghost-padded *local* extent of direction 0. The layout that `build_site`
expects, and that every sibling uses, is
`FixPosition + Orbital * prod_d Lt[d]` -- see `SimulationLMU.cpp:175-178` and
`Custom/SimulationLDOSOperator.cpp:111`. `custom_two_local` gets it right.

The check is against `ldos`, which computes the same moments by a verified route.
The lattice has to break the A/B symmetry, or the two sublattices agree trivially
and the test passes while measuring nothing: on-site +1 and -1 with
`spectrum_range=[-6, 6]` puts mu_1 = <j|Htilde|j> at +-1/6, which separates them.

Both runs use the *same* identity vertex, so the only thing that differs between
them is the requested sublattice. Without that the vertex and the seeding are
confounded and a failure cannot be attributed.

Tracked as BUG-009. The A case is asserted normally -- it passes today and would
catch a regression -- and the B case is xfail, because the fix lives in `Src/`
and is the owner's call under INV-003.
"""

import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
KITEX = REPO / "build" / "KITEx"

pytestmark = pytest.mark.skipif(not KITEX.exists(), reason="needs a built KITEx")


def _run(tmp_path, which, route):
    """Build a config for one sublattice by one route, run KITEx, return the moments."""
    import kite
    from kite import custom
    from kite import lattice as latt

    lattice = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lattice.add_sublattices(['A', [0, 0], +1.0], ['B', [0.5, 0.5], -1.0])
    for relative_index in ([0, 0], [-1, 0], [0, -1], [-1, -1]):
        lattice.add_hoppings([relative_index, 'A', 'B', -1.0])

    configuration = kite.Configuration(
        divisions=[1, 1], length=[16, 16], boundaries=["periodic"] * 2,
        is_complex=True, precision=1, spectrum_range=[-6, 6])
    calculation = kite.Calculation(configuration)
    energy = np.linspace(-0.1, 0.1, 3)

    if route == "ldos":
        calculation.ldos(energy=energy, num_moments=8, num_disorder=1,
                         position=[[0, 0]], sublattice=[which])
        dataset = "/Calculation/ldos/lMU"
    else:
        calculation.add_orbital_index('A', 0)
        calculation.add_orbital_index('B', 1)
        # Identity vertex: identical in both runs, so only the seed differs.
        calculation.add_orbital_coupling('A', 'A', 1.0, 'l0')
        calculation.add_orbital_coupling('B', 'B', 1.0, 'l0')
        calculation.custom_one_local(
            stream_=custom.Vertex(8, [[1.0, "l0"]]), energy_=energy,
            num_disorder_=1, position_=[[0, 0]], sublattice_=[which])
        dataset = "/Calculation/CustomOneLocal/Gamma"

    path = tmp_path / f"{route}_{which}.h5"
    kite.config_system(lattice, configuration, calculation, filename=str(path))
    subprocess.run([str(KITEX), path.name], cwd=str(tmp_path),
                   capture_output=True, check=True)

    import h5py
    with h5py.File(path, "r") as handle:
        return np.real(np.asarray(handle[dataset])).ravel()[:4]


def test_the_two_sublattices_are_actually_distinguishable(tmp_path):
    """Guard the guard: if this fails, the test below proves nothing."""
    a = _run(tmp_path, "A", "ldos")
    b = _run(tmp_path, "B", "ldos")
    assert a[1] == pytest.approx(+1 / 6, abs=1e-5), f"A's mu_1 is {a[1]}, not +1/6"
    assert b[1] == pytest.approx(-1 / 6, abs=1e-5), f"B's mu_1 is {b[1]}, not -1/6"


def test_custom_one_local_agrees_with_ldos_on_sublattice_a(tmp_path):
    assert _run(tmp_path, "A", "custom") == pytest.approx(
        _run(tmp_path, "A", "ldos"), abs=1e-9)


@pytest.mark.xfail(reason="BUG-009: the seed index uses r.Ld[0], the ghost-padded "
                          "local extent, so sublattice B returns A's moments",
                   strict=True)
def test_custom_one_local_agrees_with_ldos_on_sublattice_b(tmp_path):
    assert _run(tmp_path, "B", "custom") == pytest.approx(
        _run(tmp_path, "B", "ldos"), abs=1e-9)
