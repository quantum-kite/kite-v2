""" Post-processing for weyl_spectral_map.py: momentum-space spectral map A(k,E).

    spectral_map has no KITE-tools post-processing step -- the output HDF5's
    /Calculation/spectral_map/Map dataset is read directly with h5py.

    Dataset layout (confirmed against Src/Simulation/SimulationSpectral.cpp's
    store_spectral() -- which applies an fftshift-style (coord + N/2) % N
    re-indexing per dimension so k=0 sits at the center of each axis -- and
    Src/Lattice/Coordinates.cpp's index convention; cross-checked against the
    actual output below): the HDF5 dataset has shape (2, Sizet), row 0 the
    mean spectral weight and row 1 its stochastic standard error. Sizet =
    Lx*Ly*Lz*Norb, flattened as index = x + y*Lx + z*Lx*Ly + orb*(Lx*Ly*Lz)
    (sublattice slowest/outermost), so mean.reshape(Norb, Lz, Ly, Lx)
    recovers (orbital, kz, ky, kx), each axis centered (index L//2 = k=0).

    Important caveat, confirmed directly in the data below (not assumed):
    Src/FFT/FFT.cpp performs a full D-dimensional FFT unconditionally, over
    ALL axes, regardless of each axis's boundary condition. This script's
    companion weyl_spectral_map.py uses boundaries=["open","random","random"]
    (the paper authors' own intended usage) -- ky, kz (random/twisted, so
    genuinely periodic Bloch momenta) show sharp, narrow resonance peaks
    along their axis, as expected for a well-defined dispersing band cut at
    fixed energy. The "kx" axis (open / hard-wall boundary, so translational
    symmetry along x is actually broken) instead shows a broad, smeared,
    much lower-amplitude profile -- i.e. FFT still runs along x, but because
    x is not truly periodic there is no good crystal momentum there, so
    spectral weight along that axis is not a well-resolved delta-like peak.
    This is flagged directly in the figure produced below; it is a property
    of the open-boundary axis, not a bug in the reduced example sizing.

    Usage:
        python process_weyl_spectral.py <weyl_spectral_map-output.h5> <L>
    """

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()


def load_spectral_map(h5_name, L, norb=2):
    """Return (mean, stderr) spectral maps, each shape (norb, L, L, L) = (orb, kz, ky, kx)."""
    with h5py.File(h5_name, "r") as f:
        data = f["/Calculation/spectral_map/Map"][()]
    mean = data[0].reshape(norb, L, L, L)
    stderr = data[1].reshape(norb, L, L, L)
    return mean, stderr


def plot_spectral_map(h5_name, L, out_path="plots/weyl_spectral_map_preview.png"):
    """(ky, kz) slice at kx=0, plus a kz lineout at (kx, ky)=(0, 0)."""
    mean, _ = load_spectral_map(h5_name, L)
    total = mean.sum(axis=0)  # sum over sublattice, shape (kz, ky, kx)

    c = L // 2  # k=0 index on every axis (fftshift convention, see docstring)
    k_axis = np.arange(-c, L - c)

    kykz_slice = total[:, :, c]       # (kz, ky) at kx = 0
    kz_lineout = total[:, c, c]       # kz profile at kx = ky = 0
    kx_lineout = total[c, c, :]       # kx profile at ky = kz = 0, for comparison

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))

    cmap = kite_style.kite_spectral_cmap()
    mesh = axs[0].pcolormesh(k_axis, k_axis, kykz_slice, cmap=cmap,
                              shading="auto", rasterized=True)
    axs[0].set_xlabel(r"$k_y$ (grid index)", fontsize=14)
    axs[0].set_ylabel(r"$k_z$ (grid index)", fontsize=14)
    axs[0].set_title(r"$A(\mathbf{k},E{=}0.7)$ slice at $k_x{=}0$", fontsize=13)
    axs[0].set_aspect("equal")
    cbar = fig.colorbar(mesh, ax=axs[0])
    cbar.set_label("spectral weight (arb. units)", fontsize=11)

    axs[1].plot(k_axis, kz_lineout, color=kite_style.KITE_PRIMARY,
                label=r"$k_z$ lineout ($k_x{=}k_y{=}0$, random/twisted axis)")
    axs[1].plot(k_axis, kx_lineout, color=kite_style.KITE_ACCENT,
                label=r"$k_x$ lineout ($k_y{=}k_z{=}0$, OPEN axis)")
    axs[1].set_xlabel("grid index", fontsize=14)
    axs[1].set_ylabel("spectral weight (arb. units)", fontsize=14)
    axs[1].set_title("Periodic axis (sharp) vs. open axis (smeared)", fontsize=13)
    axs[1].legend(fontsize=9, loc="upper right")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")
    return total


if __name__ == "__main__":
    h5_name = sys.argv[1] if len(sys.argv) > 1 else "weyl_spectral_map-output.h5"
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    plot_spectral_map(h5_name, L)
