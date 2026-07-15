"""Real-space unit-cell and reciprocal-space Brillouin-zone plotting for a native
`kite.lattice.Lattice`.

This module is intentionally separate from `kite.lattice` (which stays importable
in headless/HDF5-generation-only contexts with no plotting dependency): it pulls in
matplotlib and is only needed when the caller actually wants a figure.

Scope of this module for now (see maintenance/native-lattice-viz-plan.md):
  - `plot_unit_cell`: 2D real-space unit cell + hoppings. 2D lattices only (exactly
    2 primitive vectors) -- no 3D support yet.
  - `plot_brillouin_zone`: 2D Brillouin-zone polygon + reciprocal vectors. 2D
    lattices only -- 1D/3D BZ plotting and k-path generation (`make_path`) are not
    implemented here yet.

Band structure / H(k) construction is explicitly out of scope for this module.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def _sublattice_colors(lattice, sublattice_colors):
    """Build (or validate) a stable name -> color mapping for a lattice's sublattices.

    Uses the qualitative 'tab10' colormap by index if the caller doesn't supply one.
    'tab10' only has 10 distinct colors; lattices with more than 10 sublattices will
    see colors repeat -- pass an explicit `sublattice_colors` dict in that case.
    """
    names = list(lattice.sublattices.keys())
    if sublattice_colors is None:
        cmap = plt.get_cmap("tab10")
        return {name: cmap(i % 10) for i, name in enumerate(names)}

    missing = [name for name in names if name not in sublattice_colors]
    if missing:
        raise ValueError(
            f"sublattice_colors is missing entries for: {sorted(missing)}"
        )
    return sublattice_colors


def plot_unit_cell(lattice, repeat=(1, 1), ax=None, sublattice_colors=None,
                    length_unit="nm"):
    """Scatter sublattice positions and draw hoppings for a 2D tight-binding lattice.

    Parameters
    ----------
    lattice : kite.lattice.Lattice
        Must have exactly 2 primitive vectors (`len(lattice.vectors) == 2`);
        3D lattices are not supported by this function yet.
    repeat : (int, int)
        How many neighboring cells to show in each of the 2 lattice-vector
        directions. `repeat=(n1, n2)` shows all cells with integer indices
        `i1 in range(-n1, n1+1)` and `i2 in range(-n2, n2+1)` -- e.g. the default
        `(1, 1)` shows the 3x3 block of cells with indices -1, 0, +1 in both
        directions, so that hoppings crossing a cell boundary are visible from
        both sides.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into; a new figure/axes is created if not given.
    sublattice_colors : dict, optional
        name -> color. Built automatically from 'tab10' (stable, up to 10
        distinct sublattices) if not supplied.
    length_unit : str, optional
        Unit string appended to the x/y axis labels (default "nm", matching the
        convention used in this repo's own graphene example). Positions are
        plotted exactly as given in `lattice.sublattices[...].position` and
        `lattice.vectors` -- pass the unit that actually matches those values.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if len(lattice.vectors) != 2:
        raise NotImplementedError(
            "plot_unit_cell only supports 2D lattices (exactly 2 primitive "
            "vectors) in this version"
        )

    a1, a2 = (np.asarray(v, dtype=float) for v in lattice.vectors)
    vectors_matrix = np.array([a1, a2])

    if ax is None:
        _, ax = plt.subplots()

    colors = _sublattice_colors(lattice, sublattice_colors)

    id_to_position = {
        sub.alias_id: np.asarray(sub.position, dtype=float)
        for sub in lattice.sublattices.values()
    }

    n1_range = range(-repeat[0], repeat[0] + 1)
    n2_range = range(-repeat[1], repeat[1] + 1)
    cell_offsets = [n1 * a1 + n2 * a2 for n1 in n1_range for n2 in n2_range]

    # Sublattice positions, one legend entry per sublattice.
    for name, sub in lattice.sublattices.items():
        pos0 = id_to_position[sub.alias_id]
        pts = np.array([pos0 + offset for offset in cell_offsets])
        ax.scatter(pts[:, 0], pts[:, 1], color=colors[name], label=name,
                   zorder=3, s=40, edgecolors="black", linewidths=0.5)

    # Hoppings: each HoppingTerms is drawn from `from_id`'s position to
    # `to_id`'s position offset by `relative_index @ vectors`, repeated at every
    # displayed cell origin so boundary-crossing bonds are visible.
    for family in lattice.hoppings.values():
        for term in family.terms:
            from_pos0 = id_to_position[term.from_id]
            to_pos0 = id_to_position[term.to_id] + np.asarray(
                term.relative_index, dtype=float
            ) @ vectors_matrix
            for offset in cell_offsets:
                p_from = from_pos0 + offset
                p_to = to_pos0 + offset
                ax.plot([p_from[0], p_to[0]], [p_from[1], p_to[1]],
                        color="0.4", linewidth=0.8, zorder=1)

    ax.set_xlabel(f"x ({length_unit})")
    ax.set_ylabel(f"y ({length_unit})")
    ax.set_aspect("equal")
    ax.legend(loc="best")
    return ax


def plot_brillouin_zone(lattice, ax=None, k_path=None, length_unit="nm"):
    """Draw the 2D Brillouin zone polygon and reciprocal lattice vectors.

    Parameters
    ----------
    lattice : kite.lattice.Lattice
        Must have exactly 2 primitive vectors; only 2D Brillouin zones can be
        plotted right now (matches `Lattice.brillouin_zone()`'s current scope).
    ax : matplotlib.axes.Axes, optional
        Axes to draw into; a new figure/axes is created if not given.
    k_path : array-like of shape (N, 2), optional
        Cartesian reciprocal-space k-points (same units as `reciprocal_vectors()`,
        i.e. 1/length_unit) to overlay as a simple connected line + markers.
        Not yet produced by any `make_path()` helper in this repo (deferred to a
        future round) -- this parameter just accepts pre-computed points.
    length_unit : str, optional
        Real-space length unit that `lattice.vectors` are given in (default
        "nm"); used only to label the k-axes as its inverse.

    Returns
    -------
    matplotlib.axes.Axes
    """
    b_vectors = lattice.reciprocal_vectors()
    if len(b_vectors) != 2:
        raise NotImplementedError(
            "plot_brillouin_zone only supports 2D lattices in this version"
        )

    vertices = lattice.brillouin_zone()  # (n_vertices, 2), CCW-ordered

    if ax is None:
        _, ax = plt.subplots()

    polygon = Polygon(vertices, closed=True, facecolor="none",
                       edgecolor="black", linewidth=1.5, zorder=2)
    ax.add_patch(polygon)

    arrow_colors = ["#1f77b4", "#d62728"]
    labels = [r"$\mathbf{b}_1$", r"$\mathbf{b}_2$"]
    for b, color, label in zip(b_vectors, arrow_colors, labels):
        ax.annotate("", xy=(b[0], b[1]), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=color, linewidth=1.5))
        ax.text(b[0] * 1.08, b[1] * 1.08, label, color=color, fontsize=11)

    all_pts = np.vstack([vertices, np.array(b_vectors), [[0.0, 0.0]]])

    if k_path is not None:
        k_path = np.asarray(k_path, dtype=float)
        ax.plot(k_path[:, 0], k_path[:, 1], color="tab:green", marker="o",
                markersize=3, linewidth=1, zorder=3, label="k-path")
        ax.legend(loc="best")
        all_pts = np.vstack([all_pts, k_path])

    ax.set_xlabel(rf"$k_x$ (1/{length_unit})")
    ax.set_ylabel(rf"$k_y$ (1/{length_unit})")
    ax.set_aspect("equal")

    span = all_pts.max(axis=0) - all_pts.min(axis=0)
    margin = 0.15 * max(span.max(), 1e-12)
    ax.set_xlim(all_pts[:, 0].min() - margin, all_pts[:, 0].max() + margin)
    ax.set_ylim(all_pts[:, 1].min() - margin, all_pts[:, 1].max() + margin)

    return ax
