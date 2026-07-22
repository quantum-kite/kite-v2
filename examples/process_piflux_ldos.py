""" Post-processing for piflux_ldos_map.py: real-space LDOS map around a vacancy.

    ldos_map has no KITE-tools post-processing step -- the output HDF5's
    /Calculation/ldos_map/Map dataset is read directly with h5py.

    Dataset layout (confirmed against Src/Simulation/SimulationLDoS.cpp's
    store_ldos() and Src/Lattice/Coordinates.cpp's index convention, and
    cross-checked against the actual output below): the HDF5 dataset has
    shape (2, Sizet), with row 0 the mean LDOS and row 1 its stochastic
    standard error (Eigen writes its Sizet-by-2 array transposed into HDF5).
    Sizet = Lx * Ly * Norb, flattened as index = x + y*Lx + orb*(Lx*Ly)
    (i.e. sublattice is the slowest-varying / outermost index), so
    mean.reshape(Norb, Ly, Lx) recovers (orbital, y, x).

    Usage:
        python process_piflux_ldos.py <piflux_ldos_map-output.h5> <Lx> <Ly>
    """

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import kite_style

kite_style.apply()


def load_ldos_map(h5_name, lx, ly, norb=2):
    """Return (mean, stderr) LDOS maps, each shape (norb, ly, lx)."""
    with h5py.File(h5_name, "r") as f:
        data = f["/Calculation/ldos_map/Map"][()]
    mean = data[0].reshape(norb, ly, lx)
    stderr = data[1].reshape(norb, ly, lx)
    return mean, stderr


def plot_vacancy_ldos(h5_name, lx, ly, half_window=24,
                       out_path="plots/piflux_ldos_map_preview.png"):
    """2D real-space LDOS map (both sublattices summed), zoomed on the vacancy."""
    mean, _ = load_ldos_map(h5_name, lx, ly)
    total = mean.sum(axis=0)  # A + B sublattice LDOS at each unit cell

    cx, cy = lx // 2, ly // 2
    x0, x1 = max(cx - half_window, 0), min(cx + half_window, lx)
    y0, y1 = max(cy - half_window, 0), min(cy + half_window, ly)
    window = total[y0:y1, x0:x1]

    # LDOS spans ~2 orders of magnitude here (background ~0.02, peak ~2.9) --
    # a linear color scale washes out everything except the brightest pixel,
    # hiding the decaying Friedel-oscillation ring structure around it. Use a
    # log color scale instead; the vacancy site itself is exactly 0 (no
    # orbital there), which LogNorm can't represent, so floor the *displayed*
    # values (not the underlying data) at a small positive value -- it then
    # just renders as the darkest color, which is the physically correct
    # takeaway (LDOS suppressed to ~0) even though the true value isn't
    # literally plottable on a log axis.
    positive = window[window > 0]
    vmin = positive.min() if positive.size else 1e-3
    vmax = window.max()
    display = np.clip(window, vmin, None)

    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = kite_style.kite_spectral_cmap()
    mesh = ax.pcolormesh(np.arange(x0, x1), np.arange(y0, y1), display,
                          cmap=cmap, shading="auto", rasterized=True,
                          norm=LogNorm(vmin=vmin, vmax=vmax))
    ax.scatter([cx], [cy], marker="x", color="white", s=80,
               linewidths=2, label="vacancy (A sublattice)")
    ax.set_xlabel("x (unit cells)", fontsize=14)
    ax.set_ylabel("y (unit cells)", fontsize=14)
    ax.set_title("Pi-flux square lattice:\nLDOS(E=0) around a single vacancy",
                 fontsize=13)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=10, framealpha=0.85)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"LDOS (arb. units, log scale)", fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")
    return total


if __name__ == "__main__":
    h5_name = sys.argv[1] if len(sys.argv) > 1 else "piflux_ldos_map-output.h5"
    lx = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    ly = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    plot_vacancy_ldos(h5_name, lx, ly)
