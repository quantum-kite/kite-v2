"""Integration test: the chiral zero mode of a single graphene vacancy.

This exercises the full pipeline -- Python export, KITEx, KITE-tools -- and
asserts the one thing about this example that is exact.

Physics
-------
On a bipartite lattice with nearest-neighbour hoppings and zero on-site
energies the Hamiltonian anticommutes with the sublattice operator, so the
spectrum is symmetric under E -> -E. Removing a single site leaves an odd
dimension, and an odd-dimensional matrix with a symmetric spectrum must have a
zero eigenvalue. That state lives entirely on the sublattice OPPOSITE the
removed atom: exact diagonalization of this lattice gives B-weight 1.000000
against 7e-30 on A.

So the acceptance criterion is the B-to-A CONTRAST at E = 0, which is exact,
and its decay with distance. Deliberately NOT asserted:

  * peak height and peak width -- the zero mode is a delta function and both
    are set by the KPM kernel, scaling with num_moments/energy_scale;
  * absolute LDOS values, which are meaningless across different num_moments
    or spectrum_range;
  * a fit of the 1/r^2 exponent -- the window spans a factor ~2.5 in r.

Sizes divisible by 3 must be avoided: the Dirac point then falls on the k-grid,
the clean lattice contributes two further zero modes, and the vacancy state
becomes an ill-defined mixture. The example uses 128.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "ldos_graphene.py"
KITEX = REPO / "build" / "KITEx"
KITETOOLS = REPO / "build" / "KITE-tools"

pytestmark = pytest.mark.skipif(
    not (KITEX.exists() and KITETOOLS.exists()),
    reason="KITEx/KITE-tools not built; run cmake && make first",
)


@pytest.fixture(scope="module")
def ldos_at_zero(tmp_path_factory):
    """Run the shipped example end to end and return the E = 0 map."""
    work = tmp_path_factory.mktemp("vacancy")
    env = {"PYTHONPATH": str(REPO / "python"), "PATH": "/usr/bin:/bin"}
    subprocess.run(["python3", str(EXAMPLE)], cwd=work, check=True, env=env,
                   capture_output=True)
    h5 = work / "graphene_lattice_ldos-output.h5"
    assert h5.exists(), "the example did not export an HDF5 file"
    subprocess.run([str(KITEX), h5.name], cwd=work, check=True, capture_output=True)
    subprocess.run([str(KITETOOLS), h5.name], cwd=work, check=True, capture_output=True)

    dat = work / "ldos0.000000.dat"
    assert dat.exists(), "KITE-tools produced no E = 0 map"
    d = np.loadtxt(dat)
    return dict(i=d[:, 0], j=d[:, 1], sub=d[:, 2].astype(int), ldos=d[:, 3])


def test_lattice_size_is_not_a_multiple_of_three():
    """3 | L puts the Dirac point on the k-grid and degenerates the zero mode."""
    src = EXAMPLE.read_text()
    assert "lx = ly = 128" in src
    assert 128 % 3 != 0


def test_zero_mode_is_polarised_on_the_opposite_sublattice(ldos_at_zero):
    """The vacancy is on A, so the E = 0 weight must sit on B."""
    d = ldos_at_zero
    a_mean = d["ldos"][d["sub"] == 0].mean()
    b_mean = d["ldos"][d["sub"] == 1].mean()
    assert b_mean / a_mean > 20, (
        f"B/A contrast is only {b_mean / a_mean:.1f}:1 -- that is a wiggle, not a "
        "sublattice-polarised zero mode. Check num_moments and spectrum_range."
    )


def test_the_asymmetry_decays_with_distance(ldos_at_zero):
    """Far from the vacancy both sublattices return to the clean value."""
    d = ldos_at_zero
    r = np.hypot(d["i"] - 64, d["j"] - 64)
    ratios = []
    for lo, hi in ((0, 1.5), (1.5, 3), (3, 5), (5, 8)):
        m = (r >= lo) & (r < hi)
        a = d["ldos"][m & (d["sub"] == 0)].mean()
        b = d["ldos"][m & (d["sub"] == 1)].mean()
        ratios.append(b / max(a, 1e-12))
    assert all(x > y for x, y in zip(ratios, ratios[1:])), (
        f"contrast must fall monotonically with distance, got {ratios}")
    assert ratios[0] > 100, "the site next to the vacancy should be strongly polarised"


def test_the_b_sublattice_excess_falls_off(ldos_at_zero):
    """|psi|^2 ~ 1/r^2, so the B-site LDOS decreases away from the vacancy."""
    d = ldos_at_zero
    r = np.hypot(d["i"] - 64, d["j"] - 64)
    B = d["sub"] == 1
    near = d["ldos"][B & (r >= 1.5) & (r < 3)].mean()
    far = d["ldos"][B & (r >= 5) & (r < 8)].mean()
    assert near > far, "the zero mode must be concentrated near the vacancy"


def test_ldos_is_non_negative(ldos_at_zero):
    """It is a density."""
    assert ldos_at_zero["ldos"].min() >= 0
