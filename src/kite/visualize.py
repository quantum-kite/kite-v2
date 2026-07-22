"""Real-space unit-cell and reciprocal-space Brillouin-zone plotting for a native
`kite.lattice.Lattice`.

This module is intentionally separate from `kite.lattice` (which stays importable
in headless/HDF5-generation-only contexts with no plotting dependency): it pulls in
matplotlib and is only needed when the caller actually wants a figure.

Scope of this module for now (see maintenance/native-lattice-viz-plan.md):
  - `plot_unit_cell`: 2D real-space unit cell + hoppings. 2D lattices only (exactly
    2 primitive vectors) -- no 3D support yet.
  - `plot_brillouin_zone`: 2D Brillouin-zone polygon + reciprocal vectors. 2D
    lattices only. Accepts an optional `k_path` (either raw points, or the
    3-tuple returned by `make_path`) to overlay.
  - `make_path`: dimension-agnostic (1D/2D/3D) piecewise-linear k-path builder
    through a sequence of high-symmetry points, uniform point density measured
    in actual Cartesian reciprocal-space distance.
  - `hamiltonian_k` / `compute_bands` / `plot_bands`: build the Bloch Hamiltonian
    directly from a lattice's own hoppings and plot the resulting tight-binding
    band structure along a k-path -- one orbital per sublattice (see
    `hamiltonian_k`'s docstring for the exact scope/gauge convention).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def make_path(*points, step=0.1, point_labels=None):
    """Build a piecewise-linear k-path through a sequence of high-symmetry points.

    Consecutive points are linearly interpolated; the number of samples drawn
    on each segment is set by that segment's actual Cartesian length divided
    by `step`, so the point *density* (points per unit reciprocal-space
    distance) is uniform across the whole path -- a segment twice as long
    (in real Cartesian k-space distance) gets twice as many points, not the
    same fixed count. This matters for anything that later estimates a band
    velocity/slope from the sampled points: non-uniform density would make
    that estimate look artificially different between segments of different
    physical length.

    Units convention (read before use): `points` must be given in the same
    Cartesian reciprocal-space coordinates returned by
    `Lattice.reciprocal_vectors()` (e.g. 1/nm if the lattice's real-space
    vectors are in nm) -- NOT fractional/reduced reciprocal-lattice
    coordinates. This function does not attempt to detect or auto-convert
    fractional input; passing fractional coordinates here will silently
    produce a k-path in the wrong units (this is the same class of bug fixed
    in commit e62c839 for the ARPES k-vector convention -- one convention,
    documented, not silently dual-accepted).

    Parameters
    ----------
    *points : array-like
        Two or more k-points (each a 1D array-like of length D, D >= 1),
        given in Cartesian reciprocal-space coordinates. Consecutive points
        must be distinct -- a zero-length segment raises `ValueError` rather
        than silently producing a divide-by-zero/NaN segment.
    step : float, optional
        Cartesian distance between consecutive sampled points (same units as
        `points`), constant along the whole path. Default 0.1. Smaller ->
        finer sampling.
    point_labels : list of str, optional
        One label per input point (`len(point_labels) == len(points)`), kept
        only as metadata for later use as band-structure x-axis tick labels
        -- this function does not itself plot anything. If omitted, the
        returned `tick_labels` is None.

    Returns
    -------
    k_points : numpy.ndarray, shape (N, D)
        Sampled points along the whole path, including every input point.
    tick_indices : list of int
        Row index into `k_points` of each input high-symmetry point, in
        input order (length == number of input points) -- for placing tick
        marks on a future band-structure plot's x-axis.
    tick_labels : list of str or None
        Echoes `point_labels` (or None if not given).
    """
    if len(points) < 2:
        raise ValueError(
            f"make_path requires at least 2 k-points, got {len(points)}"
        )
    if step <= 0:
        raise ValueError(f"step must be a positive Cartesian distance, got {step}")
    if point_labels is not None and len(point_labels) != len(points):
        raise ValueError(
            "point_labels must have exactly one label per point "
            f"({len(points)} points, {len(point_labels)} labels given)"
        )

    k_points_in = [np.atleast_1d(np.asarray(p, dtype=float)) for p in points]
    dim = k_points_in[0].shape[0]
    if any(p.ndim != 1 or p.shape[0] != dim for p in k_points_in):
        raise ValueError(
            "All k-points passed to make_path must be 1D array-likes of the "
            "same length (dimensionality)"
        )

    segments = []
    tick_indices = [0]
    for k_start, k_end in zip(k_points_in[:-1], k_points_in[1:]):
        distance = np.linalg.norm(k_end - k_start)
        if distance == 0:
            raise ValueError(
                "make_path received two identical consecutive k-points "
                f"({k_start.tolist()}); zero-length segments are not allowed"
            )
        num_steps = max(int(np.ceil(distance / step)), 1)
        segments.append(np.linspace(k_start, k_end, num_steps, endpoint=False))
        tick_indices.append(tick_indices[-1] + num_steps)

    segments.append(k_points_in[-1][np.newaxis, :])
    k_points_out = np.vstack(segments)

    return k_points_out, tick_indices, point_labels


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

    # Sublattice positions may carry more components than the 2 in-plane
    # lattice vectors (e.g. a buckled monolayer like phosphorene, where each
    # site has an out-of-plane z-offset alongside its in-plane x, y). Only the
    # first 2 components (x, y) are ever plotted here -- this is a 2D
    # top-down projection -- but an in-plane cell-translation offset (built
    # from `a1`/`a2`, which are always exactly 2-component for this 2D-only
    # function) must still be added elementwise to every position component,
    # so it's zero-padded out to the position's own dimensionality first
    # rather than raising a shape-mismatch error on any position with a z
    # component.
    position_dim = len(next(iter(id_to_position.values())))
    pad_width = position_dim - 2
    if pad_width < 0:
        raise ValueError(
            "plot_unit_cell requires sublattice positions with at least 2 "
            f"components (x, y); got {position_dim}"
        )

    def _pad(in_plane_vec):
        if pad_width == 0:
            return in_plane_vec
        return np.concatenate([in_plane_vec, np.zeros(pad_width)])

    n1_range = range(-repeat[0], repeat[0] + 1)
    n2_range = range(-repeat[1], repeat[1] + 1)
    cell_offsets = [_pad(n1 * a1 + n2 * a2) for n1 in n1_range for n2 in n2_range]

    # Sublattice positions, one legend entry per sublattice. Multiple
    # sublattices can share the exact same in-plane (x, y) position -- e.g.
    # a multi-orbital site modeled as one sublattice per orbital (TMD
    # `{metal}_0/_1/_2`, all at the metal atom's single physical position).
    # Drawing those as same-size opaque markers in dict order would let the
    # last one drawn silently hide the others behind it, even though the
    # legend claims all are present. Group by (rounded) base position and
    # draw coincident groups as nested rings (largest first, at the lowest
    # zorder within the group) instead, so every sublattice stays visible --
    # without displacing any physical position, which would misrepresent the
    # actual site geometry.
    base_xy = {
        name: tuple(np.round(id_to_position[sub.alias_id][:2], 9))
        for name, sub in lattice.sublattices.items()
    }
    coincidence_groups = {}
    for name, xy in base_xy.items():
        coincidence_groups.setdefault(xy, []).append(name)

    base_marker_size = 40
    for names_at_site in coincidence_groups.values():
        names_at_site = sorted(names_at_site)  # deterministic draw order
        n_at_site = len(names_at_site)
        for i, name in enumerate(names_at_site):
            sub = lattice.sublattices[name]
            pos0 = id_to_position[sub.alias_id]
            pts = np.array([pos0 + offset for offset in cell_offsets])
            if n_at_site > 1:
                size = max(base_marker_size - i * 14, 10)
                z = 3 + i
            else:
                size = base_marker_size
                z = 3
            ax.scatter(pts[:, 0], pts[:, 1], color=colors[name], label=name,
                       zorder=z, s=size, edgecolors="black", linewidths=0.5)

    # Hoppings: each HoppingTerms is drawn from `from_id`'s position to
    # `to_id`'s position offset by `relative_index @ vectors`, repeated at every
    # displayed cell origin so boundary-crossing bonds are visible.
    for family in lattice.hoppings.values():
        for term in family.terms:
            from_pos0 = id_to_position[term.from_id]
            to_pos0 = id_to_position[term.to_id] + _pad(np.asarray(
                term.relative_index, dtype=float
            ) @ vectors_matrix)
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
    k_path : array-like of shape (N, 2), or the 3-tuple returned by `make_path`, optional
        Cartesian reciprocal-space k-points (same units as `reciprocal_vectors()`,
        i.e. 1/length_unit) to overlay as a simple connected line + markers.
        Either pass raw points directly (shape (N, 2)), or pass
        `make_path(...)`'s full `(k_points, tick_indices, tick_labels)` return
        value directly -- in the latter case the high-symmetry points are
        additionally marked with square markers (and labeled, if
        `tick_labels` was given to `make_path`).
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
        if isinstance(k_path, tuple) and len(k_path) == 3:
            k_points_path, tick_indices, tick_labels = k_path
        else:
            k_points_path, tick_indices, tick_labels = k_path, None, None

        k_points_path = np.asarray(k_points_path, dtype=float)
        ax.plot(k_points_path[:, 0], k_points_path[:, 1], color="tab:green",
                marker="o", markersize=3, linewidth=1, zorder=3, label="k-path")

        if tick_indices is not None:
            ticks = k_points_path[tick_indices]
            ax.scatter(ticks[:, 0], ticks[:, 1], color="tab:green", marker="s",
                       s=50, zorder=4, edgecolors="black", linewidths=0.8)
            if tick_labels is not None:
                for pt, label in zip(ticks, tick_labels):
                    ax.annotate(label, xy=(pt[0], pt[1]), xytext=(3, 3),
                                textcoords="offset points", fontsize=10, zorder=5)

        ax.legend(loc="best")
        all_pts = np.vstack([all_pts, k_points_path])

    ax.set_xlabel(rf"$k_x$ (1/{length_unit})")
    ax.set_ylabel(rf"$k_y$ (1/{length_unit})")
    ax.set_aspect("equal")

    span = all_pts.max(axis=0) - all_pts.min(axis=0)
    margin = 0.15 * max(span.max(), 1e-12)
    ax.set_xlim(all_pts[:, 0].min() - margin, all_pts[:, 0].max() + margin)
    ax.set_ylim(all_pts[:, 1].min() - margin, all_pts[:, 1].max() + margin)
    return ax


def hamiltonian_k(lattice, k):
    """Build the Bloch tight-binding Hamiltonian H(k) directly from a Lattice's hoppings.

    Scope: `kite.lattice.Lattice.add_one_sublattice`/`add_one_hopping` only ever store a
    scalar onsite energy / hopping amplitude (each wrapped as a 1x1 array internally) --
    there is no multi-orbital-per-sublattice matrix support at this level (that lives in
    the separate `kite.custom` orbital-coupling machinery, used for building observables
    like spin/orbital operators, not the base tight-binding Hamiltonian). Consequently one
    sublattice == one orbital here: H(k) has shape (lattice.nsub, lattice.nsub), indexed by
    `Sublattice.alias_id`.

    Gauge / sign convention (derived, not assumed -- see maintenance/native-lattice-viz-plan.md
    section 5): matches KITE's own ARPES plane-wave state
    |k> = (1/sqrt(N)) sum_{r,alpha} w_alpha * exp(i k.(r + d_alpha)) |r,alpha>
    (docs/documentation/examples/spectral_function.md), which carries the FULL atomic
    position (unit cell + sublattice offset), not just the cell index. Fourier-transforming
    a HoppingFamily entry (relative_index=R, from_sub=a, to_sub=b, energy=t) with this same
    phase gives, with a MINUS sign in the exponent (not the naive plus -- get this backwards
    and centrosymmetric bands still look right while anything valley/Rashba/Weyl-asymmetric
    comes out mirror-flipped):

        H[b, a](k) += t * exp(-i k . (R_cartesian + d_b - d_a))

    `add_one_hopping` only ever stores one direction of a bond ("Does not add the complex
    conjugate hopping", per its own docstring) -- this function builds the one-directional
    sum H0(k) from exactly the stored hoppings, then returns H0(k) + H0(k)^dagger + onsite,
    which reproduces the missing reverse-direction terms automatically (the same effective
    result as `src/kite/__init__.py`'s own explicit reverse-hopping generation for the real
    HDF5 export, just built directly in k-space here rather than in real space per-bond).

    Parameters
    ----------
    lattice : kite.lattice.Lattice
    k : array-like, length D
        A single Cartesian reciprocal-space k-point (same units as
        `lattice.reciprocal_vectors()`'s output), D matching `lattice.vectors`' dimensionality.

    Returns
    -------
    numpy.ndarray, shape (nsub, nsub), complex
        Hermitian by construction; verified numerically before returning
        (raises `RuntimeError` rather than silently returning a broken matrix).
    """
    k = np.asarray(k, dtype=float)
    n = lattice.nsub
    H0 = np.zeros((n, n), dtype=complex)
    onsite = np.zeros((n, n), dtype=complex)

    for sub in lattice.sublattices.values():
        e = np.asarray(sub.energy)
        onsite[sub.alias_id, sub.alias_id] = e[0, 0]

    vectors_matrix = np.array(lattice.vectors, dtype=float)  # (m, D)

    for hop in lattice.hoppings.values():
        t = np.asarray(hop.energy)[0, 0]
        d_a = np.asarray(lattice.sublattices[hop.from_sub].position, dtype=float)
        d_b = np.asarray(lattice.sublattices[hop.to_sub].position, dtype=float)
        for term in hop.terms:
            R_cartesian = np.asarray(term.relative_index, dtype=float) @ vectors_matrix
            # Sublattice positions may carry more components than the lattice vectors do
            # (e.g. a buckled monolayer like phosphorene, where each site has an
            # out-of-plane z-offset with no corresponding periodic direction / reciprocal
            # k-component). R_cartesian is zero-padded to match d_b - d_a's dimensionality
            # (translating by a whole lattice vector never changes an out-of-plane offset),
            # then the dot product with k is truncated to k's own length -- k has no
            # component along a non-periodic direction, so that part of the position
            # difference simply contributes no phase (same top-down-projection convention
            # used by plot_unit_cell for buckled lattices).
            diff = d_b - d_a
            pad = len(diff) - len(R_cartesian)
            if pad > 0:
                R_cartesian = np.concatenate([R_cartesian, np.zeros(pad)])
            total = R_cartesian + diff
            phase = np.exp(-1j * np.dot(k, total[:len(k)]))
            H0[term.to_id, term.from_id] += t * phase

    H = H0 + H0.conj().T + onsite

    if not np.allclose(H, H.conj().T, atol=1e-10):
        raise RuntimeError(
            "hamiltonian_k built a non-Hermitian matrix -- this indicates a bug in the "
            "lattice's hopping definitions or in this function itself, not something to "
            "paper over by symmetrizing after the fact."
        )
    return H


def compute_bands(lattice, k_path):
    """Diagonalize `hamiltonian_k` at every point of a k-path.

    Parameters
    ----------
    lattice : kite.lattice.Lattice
    k_path : array-like, shape (N, D), or the 3-tuple returned by `make_path`
        Cartesian reciprocal-space k-points; pass either raw points or
        `make_path(...)`'s full return value directly.

    Returns
    -------
    numpy.ndarray, shape (N, lattice.nsub), real
        Eigenvalues at each k-point, ascending order along axis 1.
    """
    if isinstance(k_path, tuple) and len(k_path) == 3:
        k_points, _, _ = k_path
    else:
        k_points = k_path
    k_points = np.asarray(k_points, dtype=float)

    n = lattice.nsub
    bands = np.empty((k_points.shape[0], n), dtype=float)
    for i, k in enumerate(k_points):
        H = hamiltonian_k(lattice, k)
        bands[i] = np.linalg.eigvalsh(H)
    return bands


def plot_bands(lattice, k_path, ax=None, ylabel="Energy (eV)"):
    """Plot the tight-binding band structure along a k-path.

    Parameters
    ----------
    lattice : kite.lattice.Lattice
    k_path : array-like, shape (N, D), or the 3-tuple returned by `make_path`
        If the `make_path` 3-tuple is passed, high-symmetry points are marked with
        vertical grid lines and labeled on the x-axis (if `point_labels` were given
        to `make_path`); otherwise the x-axis is a plain 0..N-1 running index.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into; a new figure/axes is created if not given.
    ylabel : str, optional
        Energy-axis label, since the energy unit is whatever the lattice's own
        hopping/onsite values are given in (not fixed by this function).

    Returns
    -------
    matplotlib.axes.Axes
    """
    if isinstance(k_path, tuple) and len(k_path) == 3:
        k_points, tick_indices, tick_labels = k_path
    else:
        k_points, tick_indices, tick_labels = k_path, None, None

    bands = compute_bands(lattice, k_points)
    x = np.arange(bands.shape[0])

    if ax is None:
        _, ax = plt.subplots()

    for band_index in range(bands.shape[1]):
        ax.plot(x, bands[:, band_index], color="tab:blue", linewidth=1.2)

    if tick_indices is not None:
        for idx in tick_indices:
            ax.axvline(idx, color="0.7", linewidth=0.7, zorder=0)
        ax.set_xticks(tick_indices)
        if tick_labels is not None:
            ax.set_xticklabels(tick_labels)
    else:
        ax.set_xlim(x[0], x[-1])

    ax.set_ylabel(ylabel)
    ax.margins(x=0)
    return ax

    return ax
