""" Post-processing for rashba_zeeman_spin_texture.py: real-space (Sx, Sy, Sz)
    spin texture around a single vacancy in a Rashba+Zeeman honeycomb lattice,
    resolved PER SUBLATTICE (A and B separately) and plotted at true atomic
    positions so the honeycomb/hexagonal symmetry is actually visible.

    Dataset layout: /Calculation/ldos_map/Map_Operators/{l0..l5}, each an
    HDF5 dataset of shape (2, NumSites) -- row 0 mean, row 1 stderr -- with
    NumSites = lx*ly = 1024, one value per UNIT CELL (see
    register_spin_operators() in the companion script):
        l0, l1, l2 = Sx_A, Sy_A, Sz_A  (nonzero only on the Aup/Adn block)
        l3, l4, l5 = Sx_B, Sy_B, Sz_B  (nonzero only on the Bup/Bdn block)
    Reshape a row as .reshape(ly, lx) to index [y, x] on the unit-cell grid.

    Why per-sublattice (this replaces an earlier, WRONG version of this
    script/companion): a single Sx/Sy/Sz matrix spanning all 4 orbitals
    (Aup,Bup,Adn,Bdn) at once collapses the two-atom honeycomb basis into one
    number per unit cell -- structurally incapable of showing hexagonal
    symmetry, since it only resolves the underlying triangular Bravais
    lattice. Splitting Sx/Sy/Sz into an A-only and a B-only operator (as the
    companion script now does) preserves the two-atom basis, so each
    sublattice's own spin density can be placed at that atom's own real
    position and the true hexagonal arrangement becomes visible.

    Coordinates: real Cartesian atomic positions, not unit-cell indices.
    Honeycomb primitive vectors a1=[a,0], a2=[a/2, a*sqrt(3)/2] (a=0.24595 nm)
    meet at 60 degrees, not 90 -- plotting raw (i,j) indices on a square axis
    would shear the pattern by up to 30 degrees, hiding exactly the symmetry
    this figure needs to show. Sublattice A sits at i*a1 + j*a2 + posA and B
    at i*a1 + j*a2 + posB, with posA=[0,-a_cc/2], posB=[0,+a_cc/2],
    a_cc=0.142 nm (matching rashba_edelstein_graphene.py, this example's
    lattice source). Because A and B are offset from each other WITHIN each
    unit cell, a single pcolormesh on the unit-cell grid (which is what the
    first version of this script did) is no longer geometrically meaningful
    once the two sublattices are resolved separately -- each sublattice gets
    its own pcolormesh call, on its own (still-regular, just offset) oblique
    grid of atomic positions.

    Color-field rendering -- REVISED (2026-07-25): earlier revisions of this
    script rendered every panel as a `pcolormesh(..., shading="gouraud")`
    color field. That was reviewed against the actual data and rejected:

      1. Sx_A/Sy_A/Sz_A and Sx_B/Sy_B/Sz_B are each ONE STOCHASTIC KPM
         ESTIMATE PER ATOM (32x32 = 1024 sites per sublattice) -- there is no
         underlying continuum field being subsampled, so "interpolating
         between samples" has no physical justification here; it invents a
         smooth field the calculation never computed.
      2. The quantity is expected to be sharply, NON-smoothly varying right
         at the features that matter (the vacancy site, and -- in the
         zigzag-ribbon companion figures -- the edges), which is exactly
         where bilinear interpolation is most misleading.
      3. Two independent pcolormesh calls (one per sublattice, since A and B
         sit on different, mutually offset oblique grids) cannot blend
         color into each other at their shared boundary, so gouraud
         shading additionally produced a visible false seam/diamond-quilting
         artifact between the two sublattices' meshes -- worse than useless.

    This version instead plots DISCRETE MARKERS at the true atomic positions
    (circles = A, squares = B), color-mapped by value, with NO interpolation
    of any kind -- the classic real-space-probe convention (STM dI/dV maps,
    tight-binding/NEGF local-DOS figures). A thin, low-alpha bond skeleton is
    drawn first (background layer) purely to make the honeycomb connectivity
    legible; it carries no data of its own. This is design choice (b) from
    the two options considered (smooth field with lattice overlay vs.
    discrete markers as primary content) -- discrete markers were chosen
    BECAUSE the signal is atom-resolved and not smooth on the lattice scale,
    not merely as a style preference.

    Quiver arrows: intentionally NOT scaled proportionally to magnitude.
    |S_perp| spans from a stochastic noise floor of ~1e-5 to a peak of
    ~0.02 -- over 3 decades -- so no single linear length scale can make
    small-but-real signals visible without also making the peak arrow
    enormous (this is exactly what produced the "gigantic arrows" complaint
    on an earlier version, which used scale_units='xy' with a fixed
    physical-length-per-unit-spin scale). Instead, magnitude and direction
    are decoupled onto two separate visual channels, the standard approach
    for vector fields with large dynamic range: magnitude -> marker color
    (log-scale, see below), direction -> a fixed-length arrow (~0.35*a_cc,
    well under the 0.142 nm interatomic spacing at every magnitude). Arrows
    are also gated to only the statistically significant sites
    (|S_perp| > n_sigma * combined stderr, using ldos_map's own reported
    per-site stderr) so a fixed-length arrow is never drawn for a site whose
    direction is actually just sampling noise.

    Log color scale: per explicit user direction, |S_perp| is rendered on a
    log-scale color norm of its absolute value (LogNorm; all magnitudes are
    already >= 0 by construction, sqrt(Sx^2+Sy^2), so "abs" here is
    automatic). This is the right scale given the multi-decade dynamic range
    above -- a linear scale renders everything off the vacancy as flat zero.
    Sz, by contrast, is NOT log-scaled: it is a signed quantity (a real sign
    change is part of the physics here -- see the module docstring/plot
    itself, Sz_B flips from a uniform negative bulk background to slightly
    positive right next to the vacancy) with a modest (~2x) dynamic range, so
    a linear, zero-centered diverging color scale is both correct (log|Sz|
    would erase the sign that is the physically interesting feature) and
    sufficient (no multi-decade range to compress).

    Usage:
        python rashba_zeeman_spin_texture_process.py rashba_zeeman_spin_texture-output.h5
"""

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm, TwoSlopeNorm

import kite_style

kite_style.apply()

# Real-space lattice panels are drawn in true Cartesian nm coordinates on an
# oblique (60-degree) lattice; a rectilinear background grid is meaningless
# here (and was part of the earlier visual clutter), so it is switched off
# panel-by-panel below rather than left on from kite.mplstyle's default.
_ARROW_OUTLINE = [pe.withStroke(linewidth=1.3, foreground="black")]

# Honeycomb geometry, matching rashba_edelstein_graphene.py / the companion
# script's lattice (lengths in nm).
A_LATTICE = 0.24595
A_CC = 0.142
A1 = np.array([A_LATTICE, 0.0])
A2 = np.array([0.5 * A_LATTICE, 0.5 * A_LATTICE * np.sqrt(3.0)])
POS_A = np.array([0.0, -A_CC / 2])
POS_B = np.array([0.0, A_CC / 2])

OPERATOR_LABELS = {
    'l0': ('A', 'Sx'), 'l1': ('A', 'Sy'), 'l2': ('A', 'Sz'),
    'l3': ('B', 'Sx'), 'l4': ('B', 'Sy'), 'l5': ('B', 'Sz'),
}


def load_spin_maps(file_path, lx=32, ly=32):
    """Return nested dict maps[sublattice]['Sx'|'Sy'|'Sz'] = (mean, stderr),
    each array shape (ly, lx)."""
    maps = {'A': {}, 'B': {}}
    with h5py.File(file_path, "r") as f:
        grp = f["Calculation"]["ldos_map"]["Map_Operators"]
        for key, (sub, comp) in OPERATOR_LABELS.items():
            data = np.array(grp[key])  # (2, NumSites)
            mean = data[0].reshape(ly, lx)
            stderr = data[1].reshape(ly, lx)
            maps[sub][comp] = (mean, stderr)
    return maps


def sublattice_centers(lx, ly, offset):
    """Unit-cell-center-plus-offset atomic coordinates, shape (ly, lx, 2)."""
    i = np.arange(lx)
    j = np.arange(ly)
    ii, jj = np.meshgrid(i, j)
    X = ii * A1[0] + jj * A2[0] + offset[0]
    Y = ii * A1[1] + jj * A2[1] + offset[1]
    return X, Y


def window_slice(lx, ly, cx, cy, half_window):
    x0, x1 = max(cx - half_window, 0), min(cx + half_window + 1, lx)
    y0, y1 = max(cy - half_window, 0), min(cy + half_window + 1, ly)
    return slice(y0, y1), slice(x0, x1)


# The 3 A-B nearest-neighbor bonds per unit cell, in cell-index offsets --
# taken directly from rashba_edelstein_graphene.py's add_hoppings() calls
# (A at cell (i,j) bonds to B at (i,j), (i+1,j-1), (i,j-1)).
NN_BOND_OFFSETS = [(0, 0), (1, -1), (0, -1)]


def draw_bond_skeleton(ax, lx, ly, cx, cy, half_window, color="0.75", lw=0.7, zorder=1):
    """Thin lines connecting each A atom to its 3 real B neighbors, within
    the plotted window -- makes the honeycomb/hexagonal arrangement of atoms
    immediately recognizable instead of just two independent point sets."""
    for j in range(max(cy - half_window - 1, 0), min(cy + half_window + 2, ly)):
        for i in range(max(cx - half_window - 1, 0), min(cx + half_window + 2, lx)):
            a_xy = i * A1 + j * A2 + POS_A
            for di, dj in NN_BOND_OFFSETS:
                bi, bj = i + di, j + dj
                if 0 <= bi < lx and 0 <= bj < ly:
                    b_xy = bi * A1 + bj * A2 + POS_B
                    ax.plot([a_xy[0], b_xy[0]], [a_xy[1], b_xy[1]],
                            color=color, lw=lw, zorder=zorder, solid_capstyle="round")


def _style_lattice_axes(ax, xlabel="x (nm)", ylabel="y (nm)"):
    """Common cosmetic setup for real-space lattice panels: equal aspect (the
    honeycomb geometry is physically meaningful and must not be sheared/
    stretched), and the Cartesian background grid from kite.mplstyle switched
    off -- a rectilinear grid has no relationship to this oblique (60-degree)
    lattice and only adds visual noise."""
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_spin_texture(file_path, lx=32, ly=32, half_window=4, n_sigma=3.0,
                       out_path="plots/rashba_zeeman_spin_texture_preview.png"):
    maps = load_spin_maps(file_path, lx, ly)
    cx, cy = lx // 2, ly // 2  # vacancy unit cell (Aup removed there)
    vac_A_xy = cx * A1 + cy * A2 + POS_A

    XA, YA = sublattice_centers(lx, ly, POS_A)
    XB, YB = sublattice_centers(lx, ly, POS_B)

    sz_a, _ = maps['A']['Sz']
    sz_b, _ = maps['B']['Sz']

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(15.5, 5.4))
    seq_cmap = kite_style.kite_spectral_cmap()
    div_cmap = "RdBu_r"  # diverging, sign of Sz is physically meaningful here;
                          # ColorBrewer-verified colorblind-safe diverging map

    # Discrete-marker convention used consistently across every panel/figure
    # in this script (and matched by the AFM zigzag-ribbon figures): circles
    # for the A sublattice, squares for B. No interpolation anywhere -- see
    # module docstring for why a smooth color field is not physically
    # justified for this atom-resolved, stochastic-KPM quantity.
    marker_kw = dict(edgecolors="0.15", linewidths=0.3)

    # --- Panel A: Sz, full lattice, both sublattices (context) --------------
    # No marker edges here: at 1024 sites/sublattice packed across the whole
    # ribbon, edge outlines alone (not the data) would dominate the panel --
    # exactly the "busy for no informational reason" problem this redesign
    # is meant to fix. Small, edge-free dots read as a density field at this
    # scale while still being literally the discrete per-atom data.
    sz_absmax = max(abs(sz_a).max(), abs(sz_b).max())
    sz_norm_full = TwoSlopeNorm(vcenter=0.0, vmin=-sz_absmax, vmax=sz_absmax)
    ax_a.scatter(XA, YA, c=sz_a, cmap=div_cmap, norm=sz_norm_full,
                 s=6, marker="o", linewidths=0)
    sm_a = ax_a.scatter(XB, YB, c=sz_b, cmap=div_cmap, norm=sz_norm_full,
                         s=6, marker="s", linewidths=0)
    ax_a.plot(*vac_A_xy, marker="x", color="k", markersize=9,
              markeredgewidth=2, linestyle="none", zorder=5)
    ax_a.set_title(r"$S_z(\mathbf{r})$: full lattice (A + B)")
    _style_lattice_axes(ax_a)
    cbar_a = fig.colorbar(sm_a, ax=ax_a, fraction=0.046, pad=0.04)
    cbar_a.set_label(r"$S_z$ ($\hbar$)")

    # --- Panel B: Sz, zoomed on the vacancy neighborhood --------------------
    # Own (local) color limits, not panel A's global ones: the interesting
    # local feature here (Sz_B rising from a ~-0.055 bulk background to
    # +0.0019 right next to the vacancy) is a small perturbation on top of a
    # much larger global range, and would be invisible under panel A's scale.
    ys, xs = window_slice(lx, ly, cx, cy, half_window)
    sz_a_win, sz_b_win = sz_a[ys, xs], sz_b[ys, xs]
    local_absmax = max(abs(sz_a_win).max(), abs(sz_b_win).max())
    sz_norm_local = TwoSlopeNorm(vcenter=0.0, vmin=-local_absmax, vmax=local_absmax)
    draw_bond_skeleton(ax_b, lx, ly, cx, cy, half_window, color="0.75", lw=0.9, zorder=1)
    ax_b.scatter(XA[ys, xs], YA[ys, xs], c=sz_a_win, cmap=div_cmap,
                 norm=sz_norm_local, s=190, marker="o", zorder=3, **marker_kw)
    sm_b = ax_b.scatter(XB[ys, xs], YB[ys, xs], c=sz_b_win, cmap=div_cmap,
                         norm=sz_norm_local, s=190, marker="s", zorder=3, **marker_kw)
    ax_b.plot(*vac_A_xy, marker="x", color="k", markersize=13,
              markeredgewidth=2.5, linestyle="none", zorder=5)
    ax_b.set_title(r"$S_z(\mathbf{r})$: zoomed on vacancy" "\n(local color scale)")
    _style_lattice_axes(ax_b)
    cbar_b = fig.colorbar(sm_b, ax=ax_b, fraction=0.046, pad=0.04)
    cbar_b.set_label(r"$S_z$ ($\hbar$)")

    # --- Panel C: in-plane (Sx, Sy), zoomed, log|S_perp| color + direction --
    sx_a, sx_a_err = maps['A']['Sx']
    sy_a, sy_a_err = maps['A']['Sy']
    sx_b, sx_b_err = maps['B']['Sx']
    sy_b, sy_b_err = maps['B']['Sy']

    mag_a_win = np.sqrt(sx_a[ys, xs]**2 + sy_a[ys, xs]**2)
    mag_b_win = np.sqrt(sx_b[ys, xs]**2 + sy_b[ys, xs]**2)
    all_mag = np.concatenate([mag_a_win[mag_a_win > 0], mag_b_win[mag_b_win > 0]])
    vmin = max(np.percentile(all_mag, 5), 1e-4)
    vmax = all_mag.max()
    mag_norm = LogNorm(vmin=vmin, vmax=vmax)

    draw_bond_skeleton(ax_c, lx, ly, cx, cy, half_window, color="0.75", lw=0.9, zorder=1)
    ax_c.scatter(XA[ys, xs], YA[ys, xs], c=np.clip(mag_a_win, vmin, None),
                 cmap=seq_cmap, norm=mag_norm, s=190, marker="o", zorder=3, **marker_kw)
    mesh_c = ax_c.scatter(XB[ys, xs], YB[ys, xs], c=np.clip(mag_b_win, vmin, None),
                           cmap=seq_cmap, norm=mag_norm, s=190, marker="s", zorder=3, **marker_kw)

    # Fixed-length direction arrows (magnitude is already encoded in marker
    # color above), gated to statistically significant sites only -- see
    # module docstring for why arrow length is decoupled from magnitude here.
    # A black stroke outline keeps the (white) arrow visible against every
    # marker color, from the darkest to the brightest end of the colormap.
    arrow_len = 0.38 * A_CC
    for X, Y, sx, sy, sx_err, sy_err in [
        (XA[ys, xs], YA[ys, xs], sx_a[ys, xs], sy_a[ys, xs], sx_a_err[ys, xs], sy_a_err[ys, xs]),
        (XB[ys, xs], YB[ys, xs], sx_b[ys, xs], sy_b[ys, xs], sx_b_err[ys, xs], sy_b_err[ys, xs]),
    ]:
        mag = np.sqrt(sx**2 + sy**2)
        comb_err = np.sqrt(sx_err**2 + sy_err**2)
        significant = mag > n_sigma * comb_err
        # avoid divide-by-zero at the vacancy site itself (mag == 0 there)
        safe_mag = np.where(mag > 0, mag, 1.0)
        u = np.where(significant, sx / safe_mag * arrow_len, 0.0)
        v = np.where(significant, sy / safe_mag * arrow_len, 0.0)
        q = ax_c.quiver(X[significant], Y[significant], u[significant], v[significant],
                         color="white", angles="xy", scale_units="xy", scale=1,
                         width=0.016, headwidth=3.2, pivot="mid", zorder=4)
        q.set_path_effects(_ARROW_OUTLINE)

    ax_c.plot(*vac_A_xy, marker="x", color="k", markersize=13,
              markeredgewidth=2.5, linestyle="none", zorder=5)
    ax_c.set_title(r"$(S_x, S_y)$: zoomed on vacancy")
    _style_lattice_axes(ax_c)
    cbar_c = fig.colorbar(mesh_c, ax=ax_c, fraction=0.046, pad=0.04)
    cbar_c.set_label(r"$|S_\perp| = \sqrt{S_x^2+S_y^2}$ ($\hbar$, log scale)")

    fig.suptitle("Real-space spin texture around a vacancy: Rashba + Zeeman honeycomb",
                 fontsize=14, fontweight="bold", y=0.99)
    fig.text(0.5, 0.905,
              "A (circles) and B (squares) sublattices resolved separately at true "
              "atomic positions — discrete markers, no interpolation",
              ha="center", fontsize=10.5, color="0.35")
    fig.text(0.5, 0.005,
              r"Right panel: marker color $=\log_{10}|S_\perp|$, $S_\perp=\sqrt{S_x^2+S_y^2}$; "
              r"arrows show direction only (fixed length), drawn where $|S_\perp| > %g\sigma$ "
              "(stochastic KPM standard error)." % n_sigma,
              ha="center", fontsize=9.5, color="0.35")
    fig.tight_layout(rect=[0, 0.03, 1, 0.85])
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    return maps


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "rashba_zeeman_spin_texture-output.h5"
    plot_spin_texture(fname)
