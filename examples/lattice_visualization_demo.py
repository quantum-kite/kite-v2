""" Getting-started demo: build a lattice, visualize it, and plot its band structure

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping |t| = 1, length in units of lattice parameter |a| = 1
    Lattice: Simple square lattice, single orbital per unit cell
    Calculation type: none -- this script never touches KITEx/KITE-tools. It is the runnable
                       backing for docs/documentation/tb_model.md's getting-started tutorial,
                       using kite.lattice.Lattice + kite.visualize (native, pybinding-free) to
                       build a lattice, plot its unit cell / Brillouin zone / k-path, and
                       compute its tight-binding band structure directly from the lattice's
                       own hoppings -- no pb.Model/pb.Solver equivalent is needed, since
                       kite.visualize.compute_bands works straight off the primitive cell.

    Physics
    -------
    E(k) = -2t(cos(kx) + cos(ky)), t = 1: E(Gamma) = -4, E(X) = 0, E(S) = +4 -- these three
    values are printed directly by this script and are what confirm the native
    kite.lattice.Lattice + kite.visualize pipeline reproduces the same result as the
    equivalent pybinding-based tutorial this page previously used (pb.Lattice + pb.Solver).
    Last updated: 22/07/2026
"""

import numpy as np
import matplotlib.pyplot as plt

from kite import lattice as latt
from kite import visualize as viz
import kite_style

kite_style.apply()


def square_lattice(t=1.0):
    a1 = np.array([1, 0])  # [a] first lattice vector
    a2 = np.array([0, 1])  # [a] second lattice vector

    lat = latt.Lattice(a1=a1, a2=a2)

    onsite = 0
    lat.add_sublattices(
        # (name, position, onsite potential)
        ('A', [0, 0], onsite)
    )

    lat.add_hoppings(
        # (relative unit cell index, site from, site to, hopping energy)
        ([1, 0], 'A', 'A', -t),
        ([0, 1], 'A', 'A', -t)
    )
    return lat


def main():
    lat = square_lattice(t=1.0)

    ax = viz.plot_unit_cell(lat)
    ax.figure.savefig("plots/getting_started_lattice.png", dpi=150, bbox_inches="tight")
    plt.close(ax.figure)

    bz = lat.brillouin_zone()
    gamma = np.array([0.0, 0.0])
    x = (bz[1] + bz[2]) / 2
    s = bz[2]

    path = viz.make_path(gamma, x, s, gamma, step=0.02,
                          point_labels=[r"$\Gamma$", "X", "S", r"$\Gamma$"])

    ax_bz = viz.plot_brillouin_zone(lat, k_path=path)
    ax_bz.figure.savefig("plots/getting_started_brillouin.png", dpi=150, bbox_inches="tight")
    plt.close(ax_bz.figure)

    fig, ax_bands = plt.subplots(figsize=(6, 5))
    viz.plot_bands(lat, path, ax=ax_bands, ylabel="Energy (t)")
    fig.savefig("plots/getting_started_bands.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    for label, k in (("Gamma", gamma), ("X", x), ("S", s)):
        e = np.linalg.eigvalsh(viz.hamiltonian_k(lat, k))
        print(f"E({label}) = {e}")
    print("Expected: E(Gamma) = -4, E(X) = 0, E(S) = +4")


if __name__ == "__main__":
    main()
