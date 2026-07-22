""" Post-processing for altermagnet_arpes.py: fixed-energy 2D pocket maps.

    examples/process_arpes.py plots A(k,E) as a 2D (k-index-along-path, E)
    color map, which is the right tool for a 1D high-symmetry k-path. This
    example instead uses a 2D (kx, ky) k-GRID (see altermagnet_arpes.py's
    module docstring for why), so the natural plot is a single fixed-energy
    slice of A(k,E) reshaped back onto the (kx, ky) grid -- i.e. a
    constant-energy Fermi-surface-like map, one per spin, side by side, so
    the 90-degree rotation between them is directly visible. This script
    duplicates process_arpes.py's small ARPES.dat parsing block, rather
    than importing/modifying it, since process_arpes.py's own top-level
    code always makes and saves its own (unrelated) k-path figure.

    Usage:
        python process_altermagnet_arpes.py <up.dat> <down.dat> <n_kgrid> <target_energy>
    """

import sys
import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()


def parse_arpes_dat(name):
    """Return (energies, arpes) parsed from a KITE-tools --ARPES -S output file.

    arpes has shape (n_energies, n_kpoints); energies has shape (n_energies,).
    Same parsing logic as process_arpes.py's process_arpes().
    """
    with open(name, "r") as f:
        whole = f.readlines()

    cutarpes, cutenergies, cutvectors = -1, -1, -1
    count, found = 0, 0
    for line in whole:
        if "ARPES:" in line:
            cutarpes = count + 1
            found += 1
        elif "Energies:" in line:
            cutenergies = count + 1
            found += 1
        elif "k-vectors:" in line:
            cutvectors = count + 1
            found += 1
        count += 1
        if found == 3:
            break

    energies = np.array([float(x) for x in whole[cutenergies:cutarpes - 1]])
    arpes = np.array([[float(y) for y in line[:-1].split(" ") if y != ""]
                       for line in whole[cutarpes:]])
    return energies, arpes


def energy_slice(energies, arpes, target_energy):
    """Return (actual_energy, 1D intensity array over k-points) nearest target_energy."""
    idx = int(np.argmin(np.abs(energies - target_energy)))
    return energies[idx], arpes[idx, :]


def plot_pocket_comparison(up_file, down_file, n_kgrid, target_energy,
                            kmax=np.pi, out_path="altermagnet_arpes_pockets.pdf"):
    """Side-by-side spin-up / spin-down constant-energy 2D (kx,ky) maps."""
    e_up, arpes_up = parse_arpes_dat(up_file)
    e_dn, arpes_dn = parse_arpes_dat(down_file)

    E_up, slice_up = energy_slice(e_up, arpes_up, target_energy)
    E_dn, slice_dn = energy_slice(e_dn, arpes_dn, target_energy)

    # k_grid() in altermagnet_arpes.py builds the flat list as
    # [ [kx, ky] for kx in kxs for ky in kys ], i.e. kx changes slowest
    # (outer loop) -- reshape to (n_kgrid, n_kgrid) with kx as axis 0.
    grid_up = slice_up.reshape(n_kgrid, n_kgrid)
    grid_dn = slice_dn.reshape(n_kgrid, n_kgrid)

    kxs = np.linspace(-kmax, kmax, n_kgrid)
    kys = np.linspace(-kmax, kmax, n_kgrid)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
    cmap = kite_style.kite_spectral_cmap()
    vmax = max(grid_up.max(), grid_dn.max())

    for ax, grid, title, E in (
        (axs[0], grid_up, "spin up", E_up),
        (axs[1], grid_dn, "spin down", E_dn),
    ):
        mesh = ax.pcolormesh(kxs, kys, grid.T, cmap=cmap, vmin=0, vmax=vmax,
                              shading="auto", rasterized=True)
        ax.set_xlabel(r"$k_x$", fontsize=20)
        ax.set_title(f"{title}  ($E={E:.3f}\\,t$)", fontsize=20)
        ax.set_aspect("equal")
    axs[0].set_ylabel(r"$k_y$", fontsize=20)
    cbar = fig.colorbar(mesh, ax=axs, shrink=0.85)
    cbar.set_label(r"$A(\mathbf{k},E)$ (arb. units)", fontsize=18)

    fig.suptitle("Altermagnet spin-split ARPES pockets: "
                  "90-degree rotation between spin channels", fontsize=16)
    plt.savefig(out_path, dpi=150)  # format inferred from out_path's extension
    plt.close(fig)
    print(f"Saved {out_path}")
    return grid_up, grid_dn, kxs, kys


if __name__ == "__main__":
    up_file = sys.argv[1] if len(sys.argv) > 1 else "arpes_up.dat"
    down_file = sys.argv[2] if len(sys.argv) > 2 else "arpes_down.dat"
    n_kgrid = int(sys.argv[3]) if len(sys.argv) > 3 else 21
    target_energy = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    plot_pocket_comparison(up_file, down_file, n_kgrid, target_energy)
