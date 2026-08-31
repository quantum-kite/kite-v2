"""Regression test for kite.visualize.hamiltonian_k()'s row/col AND phase-sign
convention.

Two distinct bugs were found and fixed here (2026-07-26), independently:

1. Row/col: hamiltonian_k() built H0[to_id, from_id] instead of
   H0[from_id, to_id] for each stored HoppingFamily term -- backwards
   relative to what config_system() actually exports and KITEx actually
   propagates (confirmed directly against Src/Vector/KPM_Vector2D.cpp's real
   matrix-vector multiply: phi0[io] += hopping(ib,io)*phiM1[dist(ib,io)],
   i.e. row=from_id).

2. Phase sign: hamiltonian_k() used exp(-i*k.total) instead of
   exp(+i*k.total) -- backwards relative to KITE's actual C++ ARPES
   plane-wave state (confirmed directly in build_planewave(),
   Src/Vector/KPM_Vector2D.cpp/KPM_Vector3D.cpp: exp(+i*k.(r+d))). A stray
   minus sign here computes H(-k) instead of H(k).

Neither bug ever affected real KITEx simulation results -- config_system()
builds the HDF5 export independently of hamiltonian_k(). Both affected the
preview tool itself, and one confirmed downstream script that duplicated the
same (buggy) convention instead of comparing against real KITEx output (see
maintenance/2026-07-26-hopping-convention-audit.md and its follow-up).

Critically, BOTH bugs are INVISIBLE to any eigenvalue-based check: H and H^T
share eigenvalues (row/col swap), and H(k) and H(-k) share the same spectrum
at any k where -k is also sampled, which every closed-form gap/node-location
check implicitly does (symmetric k-paths). This is exactly how both bugs
went undetected through multiple "verified against known models" checks.

The two bugs also partially mask each other on a PAIRED hopping fixture
(both A->B and B->A stored explicitly, as in the Weyl model below): the
phase dependence collapses to cosines, which cannot distinguish exp(+ik)
from exp(-ik). An UNPAIRED complex hopping with distinct sublattice
positions (test 2 below) is required to catch the phase-sign bug --
this is why the original version of this test (row-only, paired fixture)
passed even before the phase-sign fix.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import kite
from kite import lattice as latt
from kite import visualize


def config_system_hamiltonian_k(lattice, k):
    """Independently reconstruct the Bloch Hamiltonian implied by
    config_system()'s own real-space hopping convention (row=from_id,
    col=to_id, +i*k.r phase, plus its auto-generated Hermitian-conjugate
    mirror terms -- see python/kite/__init__.py's config_system(), around the
    lattice.hoppings iteration), WITHOUT calling hamiltonian_k() itself.
    This is the ground-truth reference both tests below check
    hamiltonian_k() against."""
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


def weyl_lattice(t=1.0):
    """Complex, from_sub != to_sub, PAIRED hopping (both directions stored
    explicitly) -- same model as examples/weyl_lt.py. Exercises the row/col
    convention but NOT the phase sign (paired-hopping phase dependence
    collapses to cosines -- see module docstring)."""
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


def unpaired_lattice():
    """Distinct sublattice positions, ONE unpaired complex hopping (only one
    direction stored -- config_system()/hamiltonian_k() must generate the
    other via the Hermitian-conjugate mirror). This is the minimal case that
    actually distinguishes exp(+ik.total) from exp(-ik.total): with only one
    direction stored, the phase does not collapse to a cosine."""
    lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
    lat.add_sublattices(("A", [0.1, 0.2], 0.0), ("B", [0.4, 0.7], 0.0))
    lat.add_one_hopping([0, 0], "A", "B", 0.37 + 0.23j)
    return lat


def check(name, lattice, k):
    H_viz = visualize.hamiltonian_k(lattice, k)
    H_ref = config_system_hamiltonian_k(lattice, k)
    diff = np.max(np.abs(H_viz - H_ref))
    diff_T = np.max(np.abs(H_viz.T - H_ref))
    passed = diff < 1e-10
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: "
          f"hamiltonian_k vs reference = {diff:.3e} "
          f"(transposed comparison = {diff_T:.3e})")
    return passed


def main():
    results = [
        check("weyl (paired, row/col only)", weyl_lattice(t=1.0), [0.3, 0.7, 1.1]),
        check("unpaired complex hopping (row/col AND phase)",
              unpaired_lattice(), [0.3, 0.7]),
    ]
    if not all(results):
        raise SystemExit(
            "FAIL: hamiltonian_k() does not match config_system()'s "
            "row=from_id/col=to_id, exp(+i k.r) convention -- the row/col "
            "or phase-sign bug is back."
        )
    print("\nAll hamiltonian_k convention checks PASSED.")


if __name__ == "__main__":
    main()
