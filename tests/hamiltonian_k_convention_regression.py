"""Regression test for kite.visualize.hamiltonian_k()'s row/col convention.

Context: hamiltonian_k() (a pure-Python helper for quick band-structure
preview/verification plots) used to build H0[to_id, from_id] instead of
H0[from_id, to_id] for each stored HoppingFamily term -- backwards relative
to what config_system() actually exports and KITEx actually propagates
(confirmed directly against Src/Vector/KPM_Vector2D.cpp's real matrix-vector
multiply: phi0[io] += hopping(ib,io)*phiM1[dist(ib,io)], i.e. row=from_id).

This bug never affected real KITEx simulation results -- config_system()
builds the HDF5 export independently of hamiltonian_k(). It only affected
the preview tool itself, and any script that duplicated hamiltonian_k's
(buggy) convention instead of comparing against real KITEx/config_system
output (see maintenance notes for one confirmed instance of that).

Critically, this bug is INVISIBLE to any eigenvalue-based check (band
energies, gaps, node locations, Kramers degeneracies): H and H^T always
share the same eigenvalues, so a transposed-but-otherwise-correct
Hamiltonian passes every closed-form spectral check while still being the
wrong physical matrix for anything eigenvector-sensitive (ARPES spectral
weight, Berry curvature/Chern number, spin/orbital textures). This is
exactly how the bug went undetected through multiple "verified against
known models" checks. This test therefore compares matrix ELEMENTS against
an independent reconstruction of config_system()'s own real-space hopping
list, not just eigenvalues, for a complex, sublattice-asymmetric model.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import kite
from kite import lattice as latt
from kite import visualize


def weyl_lattice(t=1.0):
    """Complex, from_sub != to_sub hopping (the pattern the row/col bug needs
    to manifest) -- same model as examples/weyl_lt.py, inlined here so this
    test has no dependency on the examples/ directory."""
    lat = latt.Lattice(a1=[1, 0, 0], a2=[0, 1, 0], a3=[0, 0, 1])
    lat.add_sublattices(("A", [0, 0, 0], 0.0), ("B", [0, 0, 0], 0.0))
    lat.add_hoppings(
        ([1, 0, 0], "A", "B", 0.5 * t),
        ([1, 0, 0], "B", "A", 0.5 * t),
        ([0, 1, 0], "A", "B", 0.5 * t * 1j),
        ([0, 1, 0], "B", "A", -0.5 * t * 1j),
        ([0, 0, 1], "A", "A", 0.5 * t),
        ([0, 0, 1], "B", "B", -0.5 * t),
    )
    return lat


def config_system_hamiltonian_k(lattice, k):
    """Independently reconstruct the Bloch Hamiltonian implied by
    config_system()'s own real-space hopping convention (row=from_id,
    col=to_id, plus its auto-generated Hermitian-conjugate mirror terms --
    see src/kite/__init__.py's config_system(), around the lattice.hoppings
    iteration), WITHOUT calling hamiltonian_k() itself. This is the
    ground-truth reference this test checks hamiltonian_k() against."""
    n = lattice.nsub
    vectors_matrix = np.array(lattice.vectors, dtype=float)
    positions = {sub.alias_id: np.array(sub.position, dtype=float)
                 for sub in lattice.sublattices.values()}

    hoppings = []
    for hop in lattice.hoppings.values():
        energy = hop.energy
        for term in hop.terms:
            ridx = np.array(term.relative_index, dtype=float)
            hoppings.append((ridx, term.from_id, term.to_id, energy[0, 0]))
            if np.linalg.norm(ridx) == 0:
                hoppings.append((ridx, term.to_id, term.from_id, np.conj(energy[0, 0])))
            else:
                hoppings.append((-ridx, term.to_id, term.from_id, np.conj(energy[0, 0])))

    H = np.zeros((n, n), dtype=complex)
    for ridx, from_id, to_id, value in hoppings:
        R_cart = ridx @ vectors_matrix
        phase_arg = np.dot(k, R_cart + positions[to_id] - positions[from_id])
        H[from_id, to_id] += value * np.exp(1j * phase_arg)
    return H


def main():
    lat = weyl_lattice(t=1.0)
    k = [0.3, 0.7, 1.1]

    H_viz = visualize.hamiltonian_k(lat, k)
    H_ref = config_system_hamiltonian_k(lat, k)

    diff = np.max(np.abs(H_viz - H_ref))
    diff_T = np.max(np.abs(H_viz.T - H_ref))
    print(f"hamiltonian_k vs config_system-equivalent reference: max diff = {diff:.3e}")
    print(f"hamiltonian_k (transposed) vs reference: max diff = {diff_T:.3e}")

    if diff > 1e-10:
        raise SystemExit(
            "FAIL: hamiltonian_k() does not match config_system()'s row=from_id/"
            "col=to_id convention -- the row/col bug is back. (For reference, "
            "the transposed comparison above shows {:.3e}: if that one is ~0 "
            "instead, hamiltonian_k has reverted to the old, backwards "
            "convention.)".format(diff_T)
        )
    print("PASS: hamiltonian_k() matches config_system()'s convention exactly.")


if __name__ == "__main__":
    main()
