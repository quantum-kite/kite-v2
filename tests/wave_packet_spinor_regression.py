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
