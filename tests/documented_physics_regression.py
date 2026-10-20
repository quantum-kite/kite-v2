"""The physics each example page states, as an assertion.

Every one of these pages tells the reader what they should see — a quantized
plateau, an antisymmetric Hall curve, an exactly empty site. Those statements
are the point of the example, and nothing checked them: a change that moved a
plateau by 10% would leave the page confidently wrong.

Each assertion is against the result the page describes, produced by the route
the page names, with a tolerance set by that example's own statistics rather
than by what happens to pass today. The examples use few random vectors
deliberately — they run on a laptop — so the stochastic trace error is percent
level, and a tolerance tighter than that would flake.

These read the output of `tests/run_all_examples.py` and skip without it rather
than re-running KITEx, which takes hours.
"""

import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SWEEP = REPO / "sweep-data"
KITE_TOOLS = REPO / "build" / "KITE-tools"


def _data(example, filename):
    path = SWEEP / example / filename
    if not path.exists():
        pytest.skip(f"needs sweep-data/{example}/; run tests/run_all_examples.py first")
    return np.loadtxt(path)


def _rerun(example, flag):
    """Regenerate a .dat from stored moments, so the assertion follows the code."""
    directory = SWEEP / example
    h5 = next(directory.glob("*.h5"), None)
    if h5 is None or not KITE_TOOLS.exists():
        pytest.skip(f"needs sweep-data/{example}/ and a built KITE-tools")
    subprocess.run([str(KITE_TOOLS), h5.name, flag], cwd=str(directory),
                   capture_output=True, text=True, check=True)


def test_kane_mele_spin_hall_plateau_is_quantized():
    """docs/documentation/examples/custom_vertex_operators.md.

    The model has no Rashba term -- every hopping is spin-diagonal -- so s_z is
    exactly conserved, the two spin sectors are independent Chern insulators
    with C = +-1, and sigma^{s_z}_{xy} is quantized at -2 across the bulk gap.
    The page's window is |E| < 0.74t.

    The tolerance is 10% rather than something tighter because the measured
    value depends on the reconstruction broadening, and KITE-tools' default is
    finer than this example's moment count supports. At eta = 0.04, matched to
    M = 256, the plateau is -2.004; at KITE-tools' default eta = 0.01 it is
    -2.098, because `M >~ 10*energy_scale/eta` (docs/documentation/optimization.md)
    asks for M >~ 4000 and the example ships 256, so the reconstruction rings.
    That is a property of the parameters, not a disagreement between routes: the
    Python reconstruction gives -2.095 at the same eta, within 0.1% of
    KITE-tools. Do not tighten this without also fixing the broadening, which
    currently cannot be set from the command line.
    """
    _rerun("kane_mele_spin_hall", "--CustomTwo")
    data = _data("kane_mele_spin_hall", "custom_two.dat")
    mu, sigma = data[:, 0], data[:, 1]
    gap = np.abs(mu) < 0.74
    plateau = sigma[gap].mean()
    assert plateau == pytest.approx(-2.0, rel=0.10), (
        f"the spin Hall plateau is {plateau:.4f}, not near the quantized -2")
    flatness = sigma[gap].std() / abs(plateau)
    assert flatness < 0.05, (
        f"the plateau is {flatness:.1%} ragged; the page calls it flat")


def test_haldane_anomalous_hall_plateau_is_quantized():
    """docs/documentation/examples/haldane.md.

    The Haldane model at these parameters is a Chern insulator with C = 1, so
    sigma_xy sits at +1 in units of e^2/h across the topological gap.

    The page quotes a relative error below 0.1%, but for a run with
    `num_random=50`; the shipped example uses 10, and the stochastic trace
    error scales as 1/sqrt(N_R D). The tolerance here is set for the example as
    shipped, not for the figure.
    """
    data = _data("dos_dccond_haldane", "condDC.dat")
    mu, sigma = data[:, 0], data[:, 1]
    gap = np.abs(mu) < 0.2
    plateau = sigma[gap].mean()
    assert plateau == pytest.approx(1.0, rel=0.05), (
        f"the anomalous Hall plateau is {plateau:.4f}, not near the quantized +1")


def test_graphene_in_a_field_gives_an_antisymmetric_hall_staircase():
    """examples/dccond_graphene_magnetic_field.py.

    sigma_xy(mu) must be odd in mu: reversing the Fermi level reverses the Hall
    response. This is a symmetry of the response rather than a fitted value, so
    it holds regardless of how well the levels are resolved -- which matters
    here, because at 5/64 of a flux quantum per hexagon only the zeroth Landau
    level is resolved and the structure above the gap is Hofstadter subbands
    that this moment count does not converge. The symmetry survives all of it.
    """
    data = _data("dccond_graphene_magnetic_field", "condDC.dat")
    mu, sigma = data[:, 0], data[:, 1]
    assert np.allclose(mu, -mu[::-1], atol=1e-6), "the Fermi grid is not symmetric"
    asymmetry = np.abs(sigma + sigma[::-1]).max() / np.abs(sigma).max()
    assert asymmetry < 0.15, (
        f"sigma_xy is not odd in mu: |s(mu)+s(-mu)| reaches {asymmetry:.1%} of "
        "the peak, so the Hall response does not reverse with the Fermi level")


def test_the_pi_flux_vacancy_site_carries_no_weight():
    """docs/documentation/examples/markov_local_maps.md.

    "the LDOS is *exactly zero* at the vacancy site itself (there is no A
    orbital there to carry weight)". Exactly zero, not small -- the orbital is
    absent from the Hamiltonian, so this is structural rather than numerical.
    """
    import h5py

    directory = SWEEP / "piflux_ldos_map"
    h5 = next(directory.glob("*.h5"), None)
    if h5 is None:
        pytest.skip("needs sweep-data/piflux_ldos_map/; run tests/run_all_examples.py first")
    with h5py.File(h5, "r") as handle:
        if "/Calculation/ldos_map/Map" not in handle:
            pytest.skip("this run carries no ldos_map dataset")
        ldos = np.abs(np.asarray(handle["/Calculation/ldos_map/Map"]))

    assert ldos.size, "the LDOS map is empty"
    assert ldos.min() == 0.0, (
        f"no site carries exactly zero weight (minimum {ldos.min():.3e}); the "
        "page says the vacancy site does")
    assert ldos.max() > 0, "the LDOS map is identically zero"
