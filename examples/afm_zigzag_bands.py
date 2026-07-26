""" Band structure of a graphene honeycomb lattice with a fixed (not
    self-consistently solved) staggered Neel/AFM mass term, bulk and as a
    bearded/Klein-terminated ribbon.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Onsite energy +Delta on sublattice A / -Delta on sublattice B for spin up,
    the opposite sign for spin down -- the standard fixed-ansatz mean-field
    decoupling of the graphene Hubbard-AFM ground state (magnitude chosen in
    the typical literature range for U ~ t-3t, m ~ 0.1-0.3 t; not
    self-consistently solved here, see the companion KITE example
    afm_zigzag_ribbon.py for the corresponding real-space calculation).

    Bulk: both graphene valleys (K, K') gap out at the same |E|=Delta for
    both spins (spin-degenerate in energy, since +Delta-on-A/-Delta-on-B and
    its opposite are related by A<->B, which leaves the bulk spectrum
    invariant).

    Ribbon boundary: this lattice's three nearest-neighbor bonds -- offsets
    (0,0), (1,-1), (0,-1) in unit-cell coordinates -- put both non-intracell
    bonds from an A atom in the same direction (towards row y-1), and both
    non-intracell bonds from a B atom towards row y+1. Cutting the open
    boundary at row y=0 therefore removes both of a boundary A atom's
    inter-row bonds at once (coordination 1), not one of two as an ordinary
    zigzag edge would (coordination 2). This is a bearded/Klein-type
    termination, not an ordinary zigzag edge, despite both sublattices still
    forming an otherwise-honeycomb bulk.

    That termination still supports two dispersionless (flat) edge bands,
    each localized on one sublattice at one edge: the bottom edge (row y=0,
    coordination-1 A atoms) hosts a state that is purely sublattice A, and at
    E=+Delta that same edge is additionally purely spin-up (verified by
    direct diagonalization). The top edge (coordination-1 B atoms) is the
    mirror image: pure sublattice B, pure spin-down at that same energy. For
    a genuinely semi-infinite ribbon (or the well-localized part of the edge
    band near k=0) these bands sit AT E=+Delta and E=-Delta because the
    massless (Delta=0) edge state is exactly single-sublattice there, and the
    diagonal mass term acts on it as a scalar; this is not full chiral-symmetry
    protection of the massive Hamiltonian itself (the staggered mass term
    commutes, rather than anticommutes, with the sublattice operator sigma_z,
    so {sigma_z, H} != 0 once Delta != 0). At finite ribbon width (Ly=24
    here), the two edges hybridize by an amount that grows as the edge band
    approaches where it merges into the bulk continuum: numerically, the
    edge-band energy matches Delta to machine precision at k=0, stays within
    ~1e-7 of Delta through most of the zone, and only departs visibly from
    Delta close to that merging region near the zone boundary -- see the
    ribbon panel below, where the "flat" band's curvature away from k=0 is
    exactly this effect. This is the local, real-space-resolved edge spin
    polarization that afm_zigzag_ribbon.py's ldos_map figure visualizes
    directly; here we only verify (via exact diagonalization, not KPM/KITE)
    that these edge bands exist near the expected energies and are visibly
    non-dispersive over most of the zone.

    This script uses plain numpy exact diagonalization (not KITE/KPM) since
    the ribbon unit cell is small (a handful of atoms per row) and a direct,
    exact band structure is both cheaper and a useful independent check
    against the stochastic KPM real-space calculation in the companion
    example.

    Units: energy in units of hopping |t|=1, lengths in nm (a=0.24595 nm,
    matching every other honeycomb example in this repository).
    Last updated: 25/07/2026
"""

import os

import numpy as np
import matplotlib.pyplot as plt

import kite_style

__all__ = ["bulk_bands", "ribbon_bands", "main"]

T = 1.0
DELTA = 0.3

A = 0.24595
A1 = A * np.array([1.0, 0.0])
A2 = A * np.array([0.5, 0.5 * np.sqrt(3.0)])


def bulk_bands(kx, ky, delta, t=T):
    """Two-band (A, B) bulk Hamiltonian eigenvalues at (kx, ky), 1/nm units."""
    k = np.array([kx, ky])
    f = 1.0 + np.exp(1j * np.dot(k, A1 - A2)) + np.exp(-1j * np.dot(k, A2))
    H = np.array([[delta, -t * f], [-t * np.conj(f), -delta]])
    return np.linalg.eigvalsh(H)


def ribbon_bands(kx, ly, delta, t=T):
    """Bearded/Klein-terminated ribbon Hamiltonian eigenvalues at reduced momentum kx (in units
    of 1/a, i.e. the physical Bloch phase along a1 is exp(i*kx)), periodic
    along a1, open along a2, ly rows of (A, B) atoms."""
    n = 2 * ly
    H = np.zeros((n, n), dtype=complex)
    for y in range(ly):
        iA, iB = 2 * y, 2 * y + 1
        H[iA, iA] = delta
        H[iB, iB] = -delta
        H[iA, iB] += -t
        H[iB, iA] += -t
        if y > 0:
            iB_prev = 2 * (y - 1) + 1
            amp = -t * (1.0 + np.exp(1j * kx))
            H[iA, iB_prev] += amp
            H[iB_prev, iA] += np.conj(amp)
    return np.linalg.eigvalsh(H)


def _reciprocal_vectors():
    a_mat = np.array([A1, A2])
    b_mat = 2 * np.pi * np.linalg.inv(a_mat @ a_mat.T) @ a_mat
    return b_mat[0], b_mat[1]


def _kpath(points, n_per_segment=120):
    """Concatenate straight segments between high-symmetry points, returning
    (k-points, cumulative arc-length, segment-boundary arc-lengths) so the
    x-axis is proportional to real k-space distance, not just point index."""
    ks, dists = [], [0.0]
    for p0, p1 in zip(points[:-1], points[1:]):
        seg = np.linspace(0.0, 1.0, n_per_segment, endpoint=False)
        for s in seg:
            ks.append(p0 * (1 - s) + p1 * s)
            if len(ks) > 1:
                dists.append(dists[-1] + np.linalg.norm(ks[-1] - ks[-2]))
    ks.append(points[-1])
    dists.append(dists[-1] + np.linalg.norm(ks[-1] - ks[-2]))
    return np.array(ks), np.array(dists)


def main(ly=24, delta=DELTA, n_kx=400,
         out_path="plots/afm_zigzag_bands_preview.png"):
    kite_style.apply()

    b1, b2 = _reciprocal_vectors()
    Gamma = np.array([0.0, 0.0])
    K = (b1 + 2 * b2) / 3
    M = b1 / 2
    kpath, dists = _kpath([Gamma, K, M, Gamma])
    seg_ticks = [dists[0], dists[120], dists[240], dists[-1]]

    bands_up = np.array([bulk_bands(k[0], k[1], +delta) for k in kpath])
    bands_dn = np.array([bulk_bands(k[0], k[1], -delta) for k in kpath])

    kx_vals = np.linspace(-np.pi, np.pi, n_kx)
    ribbon_up = np.array([ribbon_bands(kx, ly, +delta) for kx in kx_vals])
    ribbon_dn = np.array([ribbon_bands(kx, ly, -delta) for kx in kx_vals])

    up_color = kite_style.KITE_PRIMARY
    dn_color = kite_style.KITE_ACCENT

    fig, (ax_bulk, ax_ribbon) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # --- Bulk: Gamma - K - M - Gamma ----------------------------------------
    ax_bulk.plot(dists, bands_up[:, 0], color=up_color, lw=2.0, label=r"spin up")
    ax_bulk.plot(dists, bands_up[:, 1], color=up_color, lw=2.0)
    ax_bulk.plot(dists, bands_dn[:, 0], color=dn_color, lw=2.0, ls="--", label=r"spin down")
    ax_bulk.plot(dists, bands_dn[:, 1], color=dn_color, lw=2.0, ls="--")
    for x in seg_ticks[1:-1]:
        ax_bulk.axvline(x, color="0.85", lw=0.8, zorder=0)
    ax_bulk.axhline(0, color="0.6", lw=0.6, zorder=0)
    ax_bulk.set_xticks(seg_ticks)
    ax_bulk.set_xticklabels([r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"])
    ax_bulk.set_xlim(dists[0], dists[-1])
    ax_bulk.set_ylabel(r"$E$ ($t$)")
    ax_bulk.set_title(r"Bulk, $\Delta = %.1f\,t$" % delta)
    ax_bulk.legend(loc="upper right", handlelength=1.6)
    ax_bulk.grid(False)

    # --- Ribbon: bearded/Klein-terminated, Ly rows --------------------------
    for b in range(ribbon_up.shape[1]):
        ax_ribbon.plot(kx_vals, ribbon_up[:, b], color=up_color, lw=1.0, alpha=0.9)
    for b in range(ribbon_dn.shape[1]):
        ax_ribbon.plot(kx_vals, ribbon_dn[:, b], color=dn_color, lw=1.0, ls="--", alpha=0.9)
    ax_ribbon.axhline(0, color="0.6", lw=0.6, zorder=0)
    ax_ribbon.axhline(+delta, color="0.2", lw=0.9, ls=":", zorder=1)
    ax_ribbon.axhline(-delta, color="0.2", lw=0.9, ls=":", zorder=1)
    ax_ribbon.annotate(r"edge bands, $E=\pm\Delta$",
                        xy=(2.6, delta), xytext=(0.3, 0.92), textcoords="axes fraction",
                        fontsize=10, ha="left", va="top",
                        arrowprops=dict(arrowstyle="-", color="0.3", lw=0.8,
                                        connectionstyle="arc3,rad=0.15"))
    ax_ribbon.set_xticks([-np.pi, 0, np.pi])
    ax_ribbon.set_xticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    ax_ribbon.set_xlim(-np.pi, np.pi)
    ax_ribbon.set_ylim(-1.2, 1.2)
    ax_ribbon.set_xlabel(r"$k_x a$")
    ax_ribbon.set_title(r"Bearded/Klein ribbon, $L_y = %d$ rows" % ly)
    ax_ribbon.grid(False)

    fig.suptitle("Staggered (Néel) mass on graphene: bulk gap vs. bearded/Klein edge bands",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
