""" Post-processing for afm_zigzag_ribbon.py: real-space Sz_A(r), Sz_B(r) map
    of a zigzag-terminated graphene ribbon with a fixed staggered Neel/AFM
    mass term, clean vs. a low (~5%) concentration of real vacancies.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Dataset layout: /Calculation/ldos_map/Map_Operators/{l0,l1}, each shape
    (2, NumSites) -- row 0 mean, row 1 stderr -- NumSites = lx*ly, one value
    per UNIT CELL: l0 = Sz_A, l1 = Sz_B (see register_sz_operators() in the
    companion script). Reshape a row as .reshape(ly, lx) to index [y, x].

    Design choice -- discrete markers, no interpolation, matching the
    (reviewed and corrected) convention now used across every real-space
    figure in this repository (see rashba_zeeman_spin_texture_process.py's
    module docstring for the full rationale): this is one stochastic KPM
    estimate per atom, expected to vary sharply (exponential edge decay,
    plus a genuine discontinuity at any vacancy site) rather than smoothly
    -- exactly the case where interpolating between samples would be
    misleading, not merely a style choice.

    Color scale -- SymLogNorm, not plain linear or plain log: Sz_A/Sz_B are
    SIGNED (a real sign difference between the two edges is the physically
    interesting result) but decay over roughly 2-3 orders of magnitude away
    from their edge, ruling out both a plain linear scale (would show only
    the outermost row of atoms and flatten everything else to a uniform
    "zero" color, hiding the decay length) and a plain log scale (erases the
    sign). SymLogNorm keeps sign, saturates at the (shared, see below) outer
    value, and still resolves the decay into the bulk.

    Shared color scale between the two panels (clean vs. vacancy) is
    mandatory here, not a default: the whole point of the comparison is
    whether the vacancy panel deviates from the clean one, which a
    per-panel-normalized color scale would silently hide.

    Candidate-vacancy markers: this particular disorder model only removes
    the A sublattice (both spin channels together, see
    afm_zigzag_ribbon.py's register/main -- add_vacancy('Aup') and
    add_vacancy('Adn') on one shared StructuralDisorder instance, so a
    "vacancy" here means the WHOLE A atom at that unit cell is gone,
    regardless of spin), and KITE's Python output does not separately record
    *which* unit cells were selected. Sites where Sz_A is exactly (to
    floating-point precision) zero are strong vacancy candidates -- but only
    where the surrounding clean-case signal is itself far from zero (i.e.
    near an edge): deep in the ribbon interior, Sz_A is ALSO ~0 simply from
    ordinary exponential edge decay, so an exact zero there is not
    attributable to a vacancy with any confidence. Candidate markers are
    therefore restricted to rows within `edge_rows` of either physical edge;
    this is a heuristic, explicitly flagged as such in the figure caption,
    not a ground-truth vacancy list.

    Usage:
        python afm_zigzag_ribbon_process.py \\
            afm_zigzag_ribbon_clean-output.h5 afm_zigzag_ribbon_vac5-output.h5
"""

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

import kite_style

kite_style.apply()

A = 0.24595
A_CC = 0.142
A1 = np.array([A, 0.0])
A2 = np.array([0.5 * A, 0.5 * A * np.sqrt(3.0)])
POS_A = np.array([0.0, -A_CC / 2])
POS_B = np.array([0.0, A_CC / 2])

DIV_CMAP = "RdBu_r"  # diverging, colorblind-safe (ColorBrewer-verified); sign
                      # of Sz is physically meaningful (opposite-edge polarization)


def load_sz_maps(file_path, lx, ly):
    with h5py.File(file_path, "r") as f:
        grp = f["Calculation"]["ldos_map"]["Map_Operators"]
        sz_a = np.array(grp["l0"])[0].reshape(ly, lx)
        sz_b = np.array(grp["l1"])[0].reshape(ly, lx)
    return sz_a, sz_b


def tiled_positions(lx, ly, offset, n_tile):
    """Real Cartesian atomic positions, periodically tiled n_tile times along
    a1 -- a faithful, non-fabricating representation of the actual periodic
    (period lx unit cells) system being simulated, not an extrapolation."""
    i = np.arange(lx * n_tile)
    j = np.arange(ly)
    ii, jj = np.meshgrid(i, j)
    X = (ii % lx) * A1[0] + jj * A2[0] + ii // lx * (lx * A1[0]) + offset[0]
    Y = (ii % lx) * A1[1] + jj * A2[1] + offset[1]
    return X, Y


def tile_values(values, n_tile):
    return np.tile(values, (1, n_tile))


def find_edge_vacancy_candidates(sz_a, ly, edge_rows=6, tol=1e-9):
    """Rows near the row-0 (bottom, A-terminated) edge only where Sz_A is
    exactly zero -- see module docstring for why this heuristic is NOT
    applied deep in the ribbon interior, and NOT near the far (row=ly-1)
    edge either: Sz_A decays to exactly 0 there from ordinary exponential
    edge decay in the CLEAN ribbon too (Sz_A lives on the A/bottom edge, not
    the B/top one), so a zero near the top edge carries no vacancy
    information at all."""
    rows, cols = np.where(np.abs(sz_a) < tol)
    mask = rows < edge_rows
    return rows[mask], cols[mask]


def plot_ribbon_comparison(clean_path, vac_path, lx=8, ly=24, n_tile=3,
                            out_path="plots/afm_zigzag_ribbon_preview.png"):
    sz_a_clean, sz_b_clean = load_sz_maps(clean_path, lx, ly)
    sz_a_vac, sz_b_vac = load_sz_maps(vac_path, lx, ly)

    XA, YA = tiled_positions(lx, ly, POS_A, n_tile)
    XB, YB = tiled_positions(lx, ly, POS_B, n_tile)

    # Shared color scale across BOTH panels -- see module docstring: the
    # comparison is only meaningful if a deviation in the vacancy panel can't
    # be hidden by a rescaled color axis.
    vmax = max(np.abs(sz_a_clean).max(), np.abs(sz_b_clean).max(),
               np.abs(sz_a_vac).max(), np.abs(sz_b_vac).max())
    linthresh = 3e-4  # below this, treat as effectively zero (bulk/noise floor)
    norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)

    fig, (ax_clean, ax_vac) = plt.subplots(1, 2, figsize=(11.5, 6.4), sharey=True)
    marker_kw = dict(cmap=DIV_CMAP, norm=norm, edgecolors="0.2", linewidths=0.25, s=42)

    for ax, sz_a, sz_b, title in [
        (ax_clean, sz_a_clean, sz_b_clean, "Clean ribbon"),
        (ax_vac, sz_a_vac, sz_b_vac, r"5% vacancy conc. (A sublattice)"),
    ]:
        ax.scatter(XA, YA, c=tile_values(sz_a, n_tile), marker="o", **marker_kw)
        sm = ax.scatter(XB, YB, c=tile_values(sz_b, n_tile), marker="s", **marker_kw)
        ax.set_title(title)
        ax.set_xlabel("x (nm)")
        ax.set_aspect("equal")
        ax.grid(False)

    # Candidate vacancy sites (heuristic, edge rows only -- see docstring).
    rows, cols = find_edge_vacancy_candidates(sz_a_vac, ly)
    # Labeled via the figure caption below, not an in-axes legend -- a legend
    # box placed inside either panel would sit on top of real lattice data at
    # this figure size, which is exactly the clutter this redesign removes
    # elsewhere.
    for t in range(n_tile):
        xs_v = (cols + t * lx) * A1[0] + rows * A2[0] + POS_A[0]
        ys_v = rows * A2[1] + POS_A[1]
        ax_vac.scatter(xs_v, ys_v, marker="x", color="k", s=70, linewidths=1.8, zorder=5)

    ax_clean.set_ylabel("y (nm)")
    cbar = fig.colorbar(sm, ax=[ax_clean, ax_vac], fraction=0.035, pad=0.03,
                         location="right")
    cbar.set_label(r"$S_z$ ($\hbar$, symmetric log scale)")

    fig.suptitle(r"Edge spin polarization in a Néel-gapped zigzag ribbon ($\Delta = 0.3\,t$)",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.925,
              "A (circles, bottom-edge sign) and B (squares, top-edge sign) sublattices at true "
              f"atomic positions; {n_tile} periods shown along $x$ (periodic boundary)",
              ha="center", fontsize=10, color="0.35")
    fig.text(0.5, 0.01,
              "Single disorder realization (no configurational averaging). Candidate vacancy "
              r"sites: exact-zero $S_z^A$ within 6 rows of the bottom (A) edge only "
              "(see script docstring for why this heuristic cannot be applied elsewhere).",
              ha="center", fontsize=8.5, color="0.4", style="italic")
    fig.tight_layout(rect=[0.0, 0.03, 0.92, 0.90])
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    clean = sys.argv[1] if len(sys.argv) > 1 else "afm_zigzag_ribbon_clean-output.h5"
    vac = sys.argv[2] if len(sys.argv) > 2 else "afm_zigzag_ribbon_vac5-output.h5"
    plot_ribbon_comparison(clean, vac)
