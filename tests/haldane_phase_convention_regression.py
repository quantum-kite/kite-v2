"""The Haldane phase convention, as the examples code it and the pages state it.

Three shipped examples build a Haldane next-nearest-neighbour term, and two of
them use the complex-conjugate assignment of the third. That is legitimate --
both are Haldane models -- but it means every page that quotes a sign has to say
which assignment it is describing, which is what T-064 is about.

The relative claim needs no numerics: with every other amplitude real, the two
lattices have complex-conjugate Hamiltonians, so their Chern numbers must be
opposite. The absolute signs do depend on conventions that have been wrong
before -- `kite.visualize.hamiltonian_k`'s `exp(+i k.r)` was a bug until
2026-07-26 per its own docstring -- and nothing pinned them, because the
existing plateau test goes through KITEx and would not see a change in
`hamiltonian_k`. This does.
"""

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _chern(lattice, grid=36):
    """Fukui-Hatsugai-Suzuki Chern number of the filled (lower) band."""
    from kite import visualize as viz

    reciprocal = np.asarray(lattice.reciprocal_vectors())[:, :2]
    states = np.empty((grid, grid, 2), complex)
    for i in range(grid):
        for j in range(grid):
            k = (i / grid) * reciprocal[0] + (j / grid) * reciprocal[1]
            _, vectors = np.linalg.eigh(viz.hamiltonian_k(lattice, k))
            states[i, j] = vectors[:, 0]

    total = 0.0
    for i in range(grid):
        for j in range(grid):
            i2, j2 = (i + 1) % grid, (j + 1) % grid
            total += np.angle(
                np.vdot(states[i, j], states[i2, j]) *
                np.vdot(states[i2, j], states[i2, j2]) *
                np.vdot(states[i2, j2], states[i, j2]) *
                np.vdot(states[i, j2], states[i, j]))
    return total / (2 * np.pi)


def test_the_reference_example_has_chern_number_plus_one():
    """`dos_dccond_haldane` codes -i*t2 on A-A; its documented plateau is +1."""
    lattice = _load("ref", "examples/dos_dccond_haldane.py").haldane()
    assert _chern(lattice) == pytest.approx(1.0, abs=1e-6), (
        "the reference Haldane example no longer has C = +1. Either the lattice "
        "changed or hamiltonian_k's exp(+i k.r) convention flipped -- the latter "
        "silently inverts every Chern sign stated in docs/")


def test_the_orbital_magnetization_example_has_the_opposite_sign():
    """It codes +i*t2 on A-A, the complex conjugate, so C = -1."""
    lattice = _load("orb", "examples/haldane_orbital_magnetization.py").haldane()
    assert _chern(lattice) == pytest.approx(-1.0, abs=1e-6), (
        "haldane_orbital_magnetization no longer has C = -1; its docstring and "
        "docs/documentation/examples/orbital_magnetization.md both state it")


def test_the_two_assignments_are_genuinely_opposite():
    """The convention-free half of the claim.

    Complex conjugation of a real-basis Hamiltonian is spinless time reversal,
    so it must flip C whatever the row/col and Fourier conventions are. If this
    ever fails, the two lattices are no longer conjugates and the pages that
    call them mirror images are wrong.
    """
    ref = _chern(_load("ref2", "examples/dos_dccond_haldane.py").haldane())
    orb = _chern(_load("orb2", "examples/haldane_orbital_magnetization.py").haldane())
    assert ref == pytest.approx(-orb, abs=1e-6), (
        f"the two Haldane examples give C = {ref:.3f} and {orb:.3f}, which are "
        "not opposite; they are supposed to be complex conjugates")


@pytest.mark.parametrize("path,must_say", [
    ("examples/kane_mele_spin_hall.py", "-i*t2 on A-A"),
    ("examples/dos_dccond_haldane.py", "-i*t2 on A-A"),
    ("examples/haldane_orbital_magnetization.py", "+i*t2 on A-A"),
    ("docs/documentation/examples/custom_vertex_operators.md", r"-i\,t_2$ on A–A"),
    ("docs/documentation/examples/haldane.md", r"$-i\,t_2$ on A–A"),
])
def test_every_page_states_the_amplitude_it_actually_codes(path, must_say):
    """T-064: one convention, stated as the amplitude rather than as a bare phi.

    A bare "phi = pi/2" is not a convention-free fact -- absorbing the leading
    minus into the exponential swaps the label -- so the pages quote the
    amplitude as coded instead.
    """
    assert must_say in (REPO / path).read_text(), (
        f"{path} no longer states the next-nearest-neighbour amplitude as "
        f"{must_say!r}; a reader cannot then tell which of the two conjugate "
        "Haldane models it builds")


def test_no_page_claims_the_two_assignments_are_equivalent():
    """They share a spectrum and differ in the sign of every chiral quantity."""
    for path in ("examples/haldane_orbital_magnetization.py",
                 "examples/kane_mele_spin_hall.py",
                 "docs/documentation/examples/haldane.md"):
        text = (REPO / path).read_text()
        assert not re.search(r"equivalent \+/- ?i|equivalent \$\\pm i", text), (
            f"{path} calls the two conjugate assignments equivalent. They have "
            "opposite Chern number, opposite sigma_xy and opposite orbital "
            "magnetization")


def test_the_fourier_sign_is_pinned_directly():
    """The Chern number cannot see this, so it needs its own test.

    Flipping `hamiltonian_k`'s `exp(+i k.r)` to `exp(-i k.r)` evaluates H(-k)
    instead of H(k). In two dimensions `k -> -k` has Jacobian +1, so it is
    orientation-PRESERVING and the Chern number is unchanged -- I measured
    C = +1 under both signs. The tests above therefore do not pin this, despite
    it being the convention that was wrong in `hamiltonian_k` until 2026-07-26.

    A single matrix element does pin it. The A-A block at wavevector k is
    `sum_R (-i t2) exp(i k.R)` over the three stored next-nearest-neighbour
    vectors, plus its conjugate from the reverse bonds; its imaginary part is
    odd under k -> -k at a generic k.
    """
    from kite import visualize as viz

    lattice = _load("fourier", "examples/dos_dccond_haldane.py").haldane()
    reciprocal = np.asarray(lattice.reciprocal_vectors())[:, :2]
    k = 0.37 * reciprocal[0] + 0.21 * reciprocal[1]      # generic, not a symmetry point

    forward = viz.hamiltonian_k(lattice, k)[0, 0]
    reversed_ = viz.hamiltonian_k(lattice, -k)[0, 0]

    assert abs(forward.imag) < 1e-12 and abs(reversed_.imag) < 1e-12, (
        "the A-A diagonal is not real; H(k) should be Hermitian")
    assert not np.isclose(forward.real, reversed_.real, atol=1e-9), (
        f"H[A,A](k) = {forward.real:.6f} equals H[A,A](-k) at a generic k, so "
        "the Haldane term is not breaking k -> -k and this test cannot pin the "
        "Fourier sign")

    # The sign itself, built from the lattice as the example defines it --
    # nothing hardcoded, because the lattice constant is 0.24595 nm, not 1.
    vectors = np.asarray(lattice.vectors)
    a1, a2 = vectors[0][:2], vectors[1][:2]
    t2 = 1.0 / 10          # the example's t2 = t/10 at its default t = 1
    expected = sum((-1j * t2) * np.exp(1j * np.dot(k, R)) +
                   np.conj(-1j * t2) * np.exp(-1j * np.dot(k, R))
                   for R in (a1, -a2, a2 - a1))
    assert forward.real == pytest.approx(expected.real, abs=1e-9), (
        f"H[A,A](k) is {forward.real:.6f} against {expected.real:.6f} built "
        "from the stored hoppings with exp(+i k.R). hamiltonian_k's Fourier "
        "sign has flipped -- it then returns H(-k), which no Chern-number test "
        "can detect in 2D")
