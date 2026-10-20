"""The sign convention custom_two applies, and the premise behind it.

`tests/custom_two_same_vertex_sign_regression.py` already checks that
`custom_two(v_x, v_x)` agrees in *sign* with `conductivity_dc(xx)`, by asserting
a correlation above 0.7. Correlation is invariant under any positive rescaling,
so it cannot see a magnitude error, and it says nothing about the vertex pairs
or disorder strengths where the sign is actually decided on a knife edge. This
file covers both.

**Why the convention needs care.** `store_custom_two` writes Tr[T_m B T_n A] --
the bare Chebyshev index on the row (Src/Simulation/Custom/SimulationRankTwo.cpp
:169) -- while the reconstruction was written against Gamma2D.cpp:138, which
puts the sandwiched index there. The two storage conventions are exact
transposes, so `customtwo.cpp` undoes the transpose when it reads Gamma. It used
instead to infer a global sign at runtime from how real Gamma was; that is the
right discriminant for a Gamma that is purely symmetric or purely antisymmetric,
but it was an L1 ratio against a 0.5 threshold measured on stochastic data, and
several shipped Gammas are neither.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

REPO = Path(__file__).resolve().parents[1]
SWEEP = REPO / "sweep-data"
KITE_TOOLS = REPO / "build" / "KITE-tools"

def _gammas():
    """Every stored custom_two Gamma in the sweep output, newest run wins."""
    found = []
    for path in sorted(SWEEP.glob("*/*.h5")):
        try:
            with h5py.File(path, "r") as handle:
                if "/Calculation/CustomTwo/Gamma" not in handle:
                    continue
                gamma = np.array(handle["/Calculation/CustomTwo/Gamma"])
        except OSError:
            continue
        found.append((path.parent.name + "/" + path.stem, gamma))
    return found


GAMMAS = _gammas()
needs_sweep = pytest.mark.skipif(
    not GAMMAS, reason="needs sweep-data/; run tests/run_all_examples.py first")


# --------------------------------------------------------------------------
# Magnitude, which the existing correlation check cannot see.
# --------------------------------------------------------------------------

GRAPHENE = (Path(__file__).resolve().parent / "fixtures" /
            "graphene_custom_two_xx_check-output.h5")


def _run_tools(flag, tmp_path):
    """Run KITE-tools for one calculation and return its (mu, sigma) columns.

    The run happens on a copy in tmp_path: KITE-tools writes its `.dat` beside
    the file it reads, and the fixture directory is tracked.
    """
    local = tmp_path / GRAPHENE.name
    if not local.exists():
        shutil.copy(GRAPHENE, local)
    subprocess.run([str(KITE_TOOLS), local.name, flag, "-F", "-1", "1", "100"],
                   cwd=str(tmp_path), capture_output=True, text=True, check=True)
    name = {"--CustomTwo": "custom_two.dat", "--CondDC": "condDC.dat"}[flag]
    data = np.loadtxt(tmp_path / name)
    return data[:, 0], data[:, 1]


@pytest.mark.skipif(not (GRAPHENE.exists() and KITE_TOOLS.exists()),
                    reason="needs the graphene fixture and a built KITE-tools")
def test_custom_two_and_conductivity_dc_agree_in_sign(tmp_path):
    """Same quantity, same file: sigma_xx must not come out negated."""
    _, custom = _run_tools("--CustomTwo", tmp_path)
    _, builtin = _run_tools("--CondDC", tmp_path)
    significant = np.abs(builtin) > 0.4 * np.abs(builtin).max()
    ratio = custom[significant] / builtin[significant]
    assert (ratio > 0).all(), (
        "custom_two(v_x, v_x) and conductivity_dc(xx) disagree in sign on the "
        "same HDF5 file -- this is the A=B case the module docstring describes")


@pytest.mark.skipif(not (GRAPHENE.exists() and KITE_TOOLS.exists()),
                    reason="needs the graphene fixture and a built KITE-tools")
def test_custom_two_is_exactly_twice_conductivity_dc(tmp_path):
    """Pins the normalization, which a correlation test cannot.

    The two paths differ by a factor of two in source: customtwo.cpp:349 uses
    `den = 2.0*num_orbitals*spin_degeneracy/...` where conductivity_dc.cpp:376
    uses `den = num_orbitals*spin_degeneracy/...`. Measured here on the same
    file, the median pointwise ratio is 1.95.

    That factor is a convention difference, not an error. The two were
    calibrated against references that differ in whether spin is in the model:
    conductivity_dc against the Haldane e^2/h plateau, which is spinless, and
    custom_two against the Kane-Mele spin Hall plateau, whose four sublattices
    (Aup, Bup, Adn, Bdn) carry spin explicitly. The lattice used here is the
    spinless one (sublattices A and B, NOrbitals = 2), so a user who wants a
    spin-degenerate result from conductivity_dc supplies the factor of 2
    themselves.

    This test pins the relationship so a change to either normalization is
    visible instead of silent.
    """
    _, custom = _run_tools("--CustomTwo", tmp_path)
    _, builtin = _run_tools("--CondDC", tmp_path)
    significant = np.abs(builtin) > 0.4 * np.abs(builtin).max()
    ratio = np.median(custom[significant] / builtin[significant])
    assert ratio == pytest.approx(2.0, rel=0.10), (
        f"custom_two/conductivity_dc median ratio is {ratio:.4f}, not ~2. One of "
        f"the two `den` constants has changed; see the docstring.")


# --------------------------------------------------------------------------
# The discriminant, and the premise underneath it.
# --------------------------------------------------------------------------

KANE_MELE = SWEEP / "kane_mele_spin_hall" / "kane_mele_spin_hall-output.h5"


@pytest.mark.skipif(not (KANE_MELE.exists() and KITE_TOOLS.exists()),
                    reason="needs the Kane-Mele fixture and a built KITE-tools")
def test_kane_mele_plateau_stays_quantized():
    """The calibration the reconstruction must not break.

    With no Rashba term the two spin sectors decouple into independent Chern
    insulators carrying +1 and -1, so sigma^{s_z}_{xy} sits on a plateau of
    magnitude 2 across the bulk gap. This is the reference the normalization was
    calibrated against, and any change to how Gamma is contracted has to
    reproduce it.
    """
    subprocess.run([str(KITE_TOOLS), KANE_MELE.name, "--CustomTwo"],
                   cwd=str(KANE_MELE.parent), capture_output=True, text=True, check=True)
    data = np.loadtxt(KANE_MELE.parent / "custom_two.dat")
    mu, sigma = data[:, 0], data[:, 1]
    plateau = sigma[np.abs(mu) < 0.6].mean()
    assert plateau == pytest.approx(-2.0, rel=0.10), (
        f"the spin Hall plateau is {plateau:.4f}, not the quantized -2")


@needs_sweep
def test_the_stored_gammas_are_not_all_cleanly_symmetric():
    """Why a single global sign cannot be inferred from the data.

    A global sign can only be correct when Gamma is purely symmetric
    (transposing is a no-op) or purely antisymmetric (transposing negates it).
    This records that the shipped Gammas are not all one or the other, which is
    why the reconstruction undoes the storage transpose directly instead of
    inferring a sign at runtime. It is a property of the stored data, so it does
    not change when the reconstruction changes -- it is here to stop the runtime
    inference being reintroduced as an apparently safe simplification.
    """
    mixed = []
    for name, gamma in GAMMAS:
        scale = np.abs(gamma).sum()
        symmetric = np.abs(gamma - gamma.T).sum() / scale
        antisymmetric = np.abs(gamma + gamma.T).sum() / scale
        if min(symmetric, antisymmetric) >= 0.25 * max(symmetric, antisymmetric):
            mixed.append(f"{name} (|G-G^T|={symmetric:.3f}, |G+G^T|={antisymmetric:.3f})")
    assert mixed, (
        "no shipped Gamma is mixed-symmetry any more; if that is genuinely true "
        "this test has lost its subject, but check the fixtures first")
