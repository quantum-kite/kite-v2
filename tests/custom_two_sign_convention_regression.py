"""The sign convention custom_two applies, and the premise behind it.

`tests/custom_two_same_vertex_sign_regression.py` already checks that
`custom_two(v_x, v_x)` agrees in *sign* with `conductivity_dc(xx)`, by asserting
a correlation above 0.7. Correlation is invariant under any positive rescaling,
so it cannot see a magnitude error, and it says nothing about the vertex pairs
or disorder strengths where the sign is actually decided on a knife edge. This
file covers both.

Some cases below fail against the current binary. Those are marked xfail
(strict) and named individually, so they report as expected failures until the
defect is fixed and then fail loudly if the markers are left behind.

**The defect.** `store_custom_two` writes Tr[T_m B T_n A] -- the bare
Chebyshev index on the row (Src/Simulation/Custom/SimulationRankTwo.cpp:169) --
while the reconstruction was written against Gamma2D.cpp:138, which puts the
sandwiched index there. The two conventions are exact transposes, so the
reconstruction needs a sign: +1 if Gamma is symmetric, -1 if antisymmetric.
Rather than undo the transpose, tools/Src/Custom/customtwo.cpp:386 infers which
case it is in by measuring how real Gamma is, against a 0.5 threshold. That
inference is sound in principle -- for a Hermitian Gamma, symmetric means real
-- but it is measured as an L1 ratio on stochastic data.
"""

import subprocess
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

REPO = Path(__file__).resolve().parents[1]
SWEEP = REPO / "sweep-data"
KITE_TOOLS = REPO / "build" / "KITE-tools"

# The L1 realness ratio the C++ compares against 0.5.
THRESHOLD = 0.5
# How far from the threshold the discriminant must sit to be trustworthy. A
# margin of 0.15 still allows a Gamma that is only 65% one way or the other.
REQUIRED_MARGIN = 0.15


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


def _realfrac(gamma):
    real = np.abs(gamma.real).sum()
    imag = np.abs(gamma.imag).sum()
    return real / (real + imag)


GAMMAS = _gammas()
needs_sweep = pytest.mark.skipif(
    not GAMMAS, reason="needs sweep-data/; run tests/run_all_examples.py first")

# The shipped cases that fail today, named individually rather than with a
# blanket xfail, so the clean cases stay genuine passes and a newly-failing one
# is a real failure rather than a silent addition to a known-bad bucket.
DISCRIMINANT_FAILS = {"kane_mele_spin_hall_disorder_W1", "kane_mele_spin_hall_disorder_W2"}
SYMMETRY_FAILS = {"kane_mele_spin_hall_disorder_W0.5",
                  "kane_mele_spin_hall_disorder_W1",
                  "kane_mele_spin_hall_disorder_W2"}


def _params(known_bad, reason):
    out = []
    for name, gamma in GAMMAS:
        stem = name.split("/")[-1].removesuffix("-output")
        marks = [pytest.mark.xfail(strict=True, reason=f"known defect: {reason}")] \
            if stem in known_bad else []
        out.append(pytest.param(name, gamma, marks=marks, id=stem))
    return out


# --------------------------------------------------------------------------
# Magnitude, which the existing correlation check cannot see.
# --------------------------------------------------------------------------

GRAPHENE = SWEEP / "graphene_custom_two_xx_check" / "graphene_custom_two_xx_check-output.h5"


def _run_tools(flag, tmp_path):
    """Run KITE-tools for one calculation and return its (mu, sigma) columns."""
    subprocess.run([str(KITE_TOOLS), GRAPHENE.name, flag, "-F", "-1", "1", "100"],
                   cwd=str(GRAPHENE.parent), capture_output=True, text=True, check=True)
    name = {"--CustomTwo": "custom_two.dat", "--CondDC": "condDC.dat"}[flag]
    data = np.loadtxt(GRAPHENE.parent / name)
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

    The two paths differ by a factor of two *in source*: customtwo.cpp:349 uses
    `den = 2.0*num_orbitals*spin_degeneracy/...` where conductivity_dc.cpp:376
    uses `den = num_orbitals*spin_degeneracy/...`, and both read the same
    `spin_degeneracy` from systemInfo. Measured here on the same file, the
    median pointwise ratio is 1.95.

    This test does NOT claim 2 is the physically correct relationship -- the two
    constants were calibrated against different references (the Kane-Mele spin
    Hall plateau and the Haldane e^2/h plateau respectively) and reconciling them
    is open. It pins the relationship so that a change to either normalization
    is visible instead of silent.
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

@needs_sweep
@pytest.mark.parametrize("name,gamma",
                         _params(DISCRIMINANT_FAILS, "discriminant sits near the threshold"))
def test_sign_discriminant_is_well_separated(name, gamma):
    """The runtime test must not be deciding a sign on a coin flip.

    customtwo.cpp:386 takes sign = +1 if realfrac > 0.5 else -1. For that to be
    trustworthy the measured value has to sit clearly on one side. On the
    shipped kane_mele_spin_hall_disorder.py it reads 0.064, 0.344, 0.438 and
    0.472 at W = 0, 0.5, 1 and 2: the antisymmetric signal decays with disorder
    while the stochastic real-part floor does not.
    """
    value = _realfrac(gamma)
    assert abs(value - THRESHOLD) > REQUIRED_MARGIN, (
        f"{name}: realfrac = {value:.4f} is within {REQUIRED_MARGIN} of the "
        f"{THRESHOLD} threshold, so the overall sign of this result is decided "
        f"by a margin comparable to the stochastic noise")


@needs_sweep
@pytest.mark.parametrize("name,gamma",
                         _params(SYMMETRY_FAILS, "Gamma is mixed-symmetry"))
def test_gamma_is_symmetric_or_antisymmetric(name, gamma):
    """The premise the whole heuristic rests on.

    A single global sign can only be correct if Gamma is purely symmetric
    (transposing changes nothing) or purely antisymmetric (transposing negates
    it). A Gamma with both parts needs them treated differently, and no choice
    of sign does that -- which is why undoing the transpose at read time is the
    exact fix and the sign branch is not.
    """
    scale = np.abs(gamma).sum()
    symmetric = np.abs(gamma - gamma.T).sum() / scale
    antisymmetric = np.abs(gamma + gamma.T).sum() / scale
    assert min(symmetric, antisymmetric) < 0.25 * max(symmetric, antisymmetric), (
        f"{name}: |G-G^T|/|G| = {symmetric:.3f} and |G+G^T|/|G| = "
        f"{antisymmetric:.3f} are comparable, so this Gamma is neither "
        f"symmetric nor antisymmetric and no single sign is correct for it")
