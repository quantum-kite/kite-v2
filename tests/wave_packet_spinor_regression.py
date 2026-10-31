"""The seeded wave packet lives in the orbital basis, and the docs say so.

`gaussian_wave_packet(spinor=...)` takes an amplitude per orbital of the unit
cell. `KPM_Vector<T,2>::build_wave_packet` writes `psi0(io, ik)` straight onto
orbital `io` under a Gaussian envelope and a plane-wave phase -- there is no
diagonalisation of H(k) anywhere in the path. So the spinor a user writes is a
statement about orbitals, not about bands, and unless it is an eigenvector of
H(k) the packet starts with weight in every band.

What is asserted here is the arithmetic a reader would do on a honeycomb
lattice, where the two orbitals are the A and B sublattices, and that the three
places documenting `spinor` say which basis it is in.
"""

import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]

#: Every place a user meets the parameter, and how to find the text that
#: describes it in each. Anchoring on the description itself matters: a whole
#: file trivially contains the words "orbital" and "band" somewhere.
def _api_table_row(text):
    for line in text.splitlines():
        if line.lstrip().startswith("| `#!python spinor`"):
            return line
    return ""


def _python_docstring(text):
    start = text.index("        spinor : np.array")
    rest = text[start + len("        spinor : np.array"):]
    end = rest.index("        width : float")
    return rest[:end]


def _page_paragraph(text):
    start = text.index("The `#!python spinor=[[1.0, 0.0]]`")
    return text[start:start + 1600]


DOCUMENTED = [
    (Path("docs/documentation/examples/time_evolution.md"), _page_paragraph),
    (Path("docs/api/kite.md"), _api_table_row),
    (Path("python/kite/__init__.py"), _python_docstring),
]


def _band_weights(hamiltonian, spinor):
    """|<band|spinor>|^2, the weight the seeded state puts in each band."""
    _, vectors = np.linalg.eigh(hamiltonian)
    return np.abs(vectors.conj().T @ np.asarray(spinor, dtype=complex)) ** 2


def _graphene_bloch(k, t=1.0):
    """The honeycomb Bloch matrix, with the hopping list the examples use."""
    deltas = [np.array([0, 0]), np.array([1, -1]), np.array([0, -1])]
    f = -t * sum(np.exp(1j * np.dot(k, d)) for d in deltas)
    return np.array([[0.0, f], [np.conj(f), 0.0]])


def _precession_bloch(k, t=1.0, field=0.1):
    """The shipped example: identical hopping on both orbitals, plus B*sigma_x.

    The on-site term couples two *different* sublattices at the same position,
    so it is off-diagonal whichever direction convention add_hoppings uses.
    """
    eps = -2 * t * (np.cos(k[0]) + np.cos(k[1]))
    return eps * np.eye(2) + field * np.array([[0.0, 1.0], [1.0, 0.0]])


def test_the_shipped_example_seeds_both_bands_equally():
    """Its precession is the beat between them, not motion within one band."""
    field = 0.1
    for k in ([0.0, 0.0], [0.7, 1.2], [-2.0, 0.9]):
        hamiltonian = _precession_bloch(np.array(k), field=field)
        values, _ = np.linalg.eigh(hamiltonian)
        assert values[1] - values[0] == pytest.approx(2 * field), (
            f"at k={k} the splitting is {values[1] - values[0]}, not 2B; the "
            "precession frequency would not be 2B")
        weights = _band_weights(hamiltonian, [1, 0])
        assert np.allclose(weights, 0.5, atol=1e-12), (
            f"at k={k} the seeded (1,0) gives band weights {weights}, not the equal "
            "split that makes the beat")


def test_the_eigenvector_spinor_selects_one_band():
    """What the page tells a reader to pass to stay in a single band."""
    for k in ([0.0, 0.0], [0.7, 1.2]):
        weights = _band_weights(_precession_bloch(np.array(k)), [1 / np.sqrt(2)] * 2)
        assert np.isclose(weights.max(), 1.0, atol=1e-12), (
            f"at k={k} (1,1)/sqrt(2) is no longer an eigenvector ({weights}); the "
            "page says it leaves <Sx> at +1/2 and nothing precesses")


def test_an_orbital_spinor_spans_every_band_on_a_honeycomb_lattice():
    """A != B sublattices: H(k) is off-diagonal, so (1,0) is never an eigenvector."""
    for k in ([0.0, 0.0], [0.5, 0.3], [2 * np.pi / 3, 2 * np.pi / 3], [1.7, -0.4]):
        weights = _band_weights(_graphene_bloch(np.array(k)), [1, 0])
        assert np.allclose(weights, 0.5, atol=1e-12), (
            f"at k={k} the orbital spinor (1,0) gives band weights {weights}, not the "
            "equal split the documentation describes")


def test_every_place_that_documents_the_spinor_names_the_basis():
    """Three doors to the same parameter; a reader may come through any of them."""
    for relative, extract in DOCUMENTED:
        path = REPO / relative
        described = extract(path.read_text())
        assert described.strip(), (
            f"{relative} no longer has a description of `spinor` where this test "
            "looks for one")
        assert re.search(r"orbital", described, re.I), (
            f"{relative} documents `spinor` without saying it is an amplitude per "
            f"orbital. What it says: {described.strip()[:200]!r}")
        assert re.search(r"\bband", described, re.I), (
            f"{relative} documents `spinor` without saying what it implies for the "
            f"bands the packet occupies. What it says: {described.strip()[:200]!r}")


# --- AC-WAVE-05: the k convention, and the spinor across the packet's spread ---

def _honeycomb():
    """The lattice the k-convention statement uses as its example."""
    from kite import lattice as latt

    a = 0.246
    lattice = latt.Lattice(a1=[a, 0], a2=[a / 2, a * np.sqrt(3) / 2])
    lattice.add_sublattices(['A', [0, 0], 0.0], ['B', [0, a / np.sqrt(3)], 0.0])
    for relative_index in ([0, 0], [-1, 0], [0, -1]):
        lattice.add_hoppings([relative_index, 'A', 'B', -1.0])
    return lattice


def test_the_documented_k_conversion_round_trips():
    """The snippet on the page is the inverse of `reciprocal_vectors()`.

    If these ever stop being inverses, the page is telling users to seed the
    wrong momentum.
    """
    lattice = _honeycomb()
    vectors = np.asarray(lattice.vectors)
    vectors = vectors[:, :vectors.shape[0]]
    reciprocal = np.asarray(lattice.reciprocal_vectors())[:, :2]

    k_frac = np.array([1 / 3, -1 / 3])
    k_cart = k_frac @ reciprocal
    assert np.allclose(k_cart @ vectors.T / (2 * np.pi), k_frac), (
        "the conversion documented in time_evolution.md does not invert "
        "reciprocal_vectors()")


def test_the_two_k_conventions_really_do_differ():
    """`hamiltonian_k` is Cartesian; `gaussian_wave_packet` is fractional.

    Asserted on eigenvalues only because that is enough to show the two are
    different states -- [1/3, -1/3] is graphene's K point, where the gap closes,
    and the same numbers read as Cartesian are nowhere near it. A test that
    passed either way would be no test at all.
    """
    from kite import visualize as viz

    lattice = _honeycomb()
    reciprocal = np.asarray(lattice.reciprocal_vectors())[:, :2]
    k_frac = np.array([1 / 3, -1 / 3])

    at_k = np.linalg.eigvalsh(viz.hamiltonian_k(lattice, k_frac @ reciprocal))
    assert np.allclose(at_k, 0, atol=1e-9), (
        f"[1/3, -1/3] converted to Cartesian is not the Dirac point: {at_k}")

    raw = np.linalg.eigvalsh(viz.hamiltonian_k(lattice, k_frac))
    assert np.abs(raw).min() > 1.0, (
        "passing the fractional pair to hamiltonian_k lands at the same place "
        "as converting it, so the page's warning describes nothing")


def test_the_spinor_is_one_row_per_k_vector():
    """The page offers per-k spinors as the remedy; the API must accept them.

    `SimulationGaussianWavePacket.cpp` sizes the spinor buffer from the row
    count of the k_vector dataset, so the two must agree.
    """
    source = (REPO / "Src" / "Simulation" / "SimulationGaussianWavePacket.cpp").read_text()
    assert re.search(r"spinor\s*=\s*Eigen::Matrix<\s*T,-1,\s*-1>::Zero\(r\.Orb,\s*dim\[0\]\)",
                     source), (
        "the C++ no longer sizes the spinor as (orbitals, number of k-vectors); "
        "the per-k spinor documented on the page may no longer be read correctly")


@pytest.mark.parametrize("must_say", [
    "fractional",
    "[1/3, -1/3]",
    "visualize-hamiltonian_k",
    "one row per $k$-vector",
])
def test_the_page_states_the_k_convention_and_the_spread(must_say):
    """AC-WAVE-05's two remaining clauses, anchored on the admonitions.

    Scoped to the two blocks that carry them rather than the whole page, so a
    passing assertion means the statement is where a reader choosing a spinor
    will meet it.
    """
    page = (REPO / "docs" / "documentation" / "examples" / "time_evolution.md").read_text()
    start = page.index('!!! warning "The spinor is held fixed')
    end = page.index("### Validation", start)
    assert must_say in page[start:end], (
        f"time_evolution.md's wave-packet admonitions no longer say {must_say!r}")


# --- the timestep convention the example documents ---------------------------

def _example_source():
    return (REPO / "examples" / "spin_precession_simple.py").read_text()


def test_the_documented_time_span_matches_the_example_it_describes():
    """The docstring quotes a span in physical time; derive it and check.

    The example asks for its step in physical units, so the span is just
    `num_points * timestep_physical` and no longer moves with `spectrum_range`.
    If anyone edits either default, the sentence silently becomes wrong --
    exactly the kind of drift a reader cannot detect.
    """
    source = _example_source()
    signature = source[source.index("def main("):source.index("output_file=")]

    step = float(re.search(r"timestep_physical=([\d.]+)", signature).group(1))
    num_points = int(re.search(r"num_points=(\d+)", signature).group(1))
    field = float(re.search(r"B=([\d.]+)", signature).group(1))

    last = round((num_points - 1) * step, 10)
    assert f"the {num_points} points span t = 0 to {last}" in source, (
        f"the docstring's arithmetic no longer matches the code: "
        f"{num_points} points of {step} span t = 0 to {last}")
    assert f"2B = {round(2 * field, 10)}" in source, (
        f"the docstring quotes the wrong precession frequency for B={field}")


def test_the_example_says_which_step_convention_it_is_using():
    """The sentence a user must read before setting either step."""
    source = _example_source()
    assert "t_rescaled = t * EnergyScale" in source, (
        "spin_precession_simple.py no longer states the relation between the "
        "rescaled step KITEx propagates in and physical time; a user who reads "
        "`timestep` as physical time gets an axis wrong by EnergyScale")
    assert "exactly one of them is given" in source, (
        "the docstring no longer says that `timestep` and `timestep_physical` "
        "are alternatives rather than two independent settings")


# --- the energy shift in the propagator --------------------------------------

def test_the_propagator_carries_the_shift_phase():
    """exp(-iHt) = exp(-i*eps0*t) * exp(-i*Htilde*ttilde); both factors must be there.

    The Bessel series supplies the second. The first is a global phase, so it
    cancels in every expectation value the driver reports -- which is why its
    absence went unnoticed. It is not unobservable in principle: it is the whole
    difference between propagating H and propagating Htilde, and anything that
    uses the amplitude rather than a bilinear in it would see it.

    Asserted structurally rather than by a run, because the behavioural
    consequence is confined to `Id`, which AC-WAVE-06 removes from the user's
    view. The phase was verified numerically when it was added: with
    spectrum_range=[-4,6] (shift 1, physical step 0.1), Id picked up exactly
    exp(-i*n*0.1) to 2e-12 while <Sz> moved by 1.5e-14 and |Id| by 6e-16.
    """
    source = (REPO / "Src" / "Simulation" / "SimulationGaussianWavePacket.cpp").read_text()
    assert '"/EnergyShift"' in source, (
        "SimulationGaussianWavePacket no longer reads /EnergyShift, so it "
        "propagates exp(-i*Htilde*ttilde) rather than exp(-i*H*t)")
    assert re.search(r"assign_value\(std::cos\(theta\),\s*-std::sin\(theta\)\)", source), (
        "the shift phase exp(-i*eps0*t) is no longer applied to the Bessel "
        "coefficients")


#: Built as a string so the two runs differ only in `spectrum_range`.
BUILD_TEMPLATE = """
import sys
sys.path.insert(0, %(examples)r)
import kite
from spin_precession_simple import simple_square_lattice_with_field
lo, hi, out = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
lattice = simple_square_lattice_with_field(t=1.0, B=0.1)
configuration = kite.Configuration(
    divisions=[1, 1], length=[64, 64], boundaries=["periodic"] * 2,
    is_complex=True, precision=1, spectrum_range=[lo, hi])
calculation = kite.Calculation(configuration)
calculation.gaussian_wave_packet(
    num_points=120, num_moments=256, num_disorder=1,
    k_vector=[[0.0, 0.0]], spinor=[[1.0, 0.0]], width=4.0, timestep=0.5,
    mean_value=[32, 32])
kite.config_system(lattice, configuration, calculation, filename=out)
"""

def test_the_shift_phase_has_the_right_sign(tmp_path):
    """Two spectrum_ranges of equal width must give the SAME propagated amplitude.

    This is the test that catches a sign error, and the reason the obvious one
    does not. Expectation values are blind to a global phase whatever its sign,
    so asserting `<Sz>` is shift-independent proves nothing about it -- that held
    just as well before the phase existed at all.

    `Id` is a bare complex amplitude, so it does see the phase. The physical
    Hamiltonian is identical in both runs (only the rescaling differs), so once
    `exp(-i*eps0*t)` is applied with the correct sign the two must agree. With
    the sign flipped the residual is ~0.25 against a signal of ~0.14 -- larger
    than the quantity itself.
    """
    import subprocess
    import sys

    import h5py

    kitex = REPO / "build" / "KITEx"
    if not kitex.exists():
        pytest.skip("needs a built KITEx")

    build = BUILD_TEMPLATE % {"examples": str(REPO / "examples")}

    amplitudes = {}
    for lo, hi, name in ((-5, 5, "sym.h5"), (-4, 6, "asym.h5")):
        subprocess.run([sys.executable, "-c", build, str(lo), str(hi), name],
                       cwd=str(tmp_path), capture_output=True, check=True)
        subprocess.run([str(kitex), name], cwd=str(tmp_path),
                       capture_output=True, check=True)
        with h5py.File(tmp_path / name, "r") as handle:
            amplitudes[name] = np.asarray(
                handle["/Calculation/gaussian_wave_packet/Id"]).ravel()

    scale = np.abs(amplitudes["sym.h5"]).max()
    residual = np.abs(amplitudes["asym.h5"] - amplitudes["sym.h5"]).max()
    assert residual < 1e-9 * scale, (
        "the propagated amplitude differs by %.3e between spectrum_range=[-5,5] "
        "and [-4,6], against a signal of %.3e. Both describe the same "
        "Hamiltonian, so exp(-i*eps0*t) is missing, has the wrong sign, or is "
        "scaled wrongly" % (residual, scale))
