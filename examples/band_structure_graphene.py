""" Tight-binding band structure of monolayer graphene, straight from the lattice definition

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV, Length in nm
    Lattice: Honeycomb (same lattice as dos_graphene.py)
    Calculation type: none -- this script never touches KITEx/KITE-tools at all. It uses
                       kite.visualize.hamiltonian_k/compute_bands/plot_bands to diagonalize
                       the Bloch Hamiltonian built directly from the lattice's own hoppings,
                       a fast (pure numpy) sanity check of a lattice definition BEFORE paying
                       for an expensive KPM/Chebyshev run -- not a replacement for KITE's own
                       stochastic machinery, just a way to catch a wrong sign/hopping/lattice
                       vector before spending compute on it.

    Physics
    -------
    Reproduces the textbook graphene tight-binding bands along the Gamma-K-M-Gamma path:
    E(Gamma) = +/-3t (all 3 NN bonds add in phase), E(M) = +/-t (van Hove saddle point),
    E(K) = 0 (the Dirac point, exactly two-fold degenerate). These three values were used
    to independently verify kite.visualize.hamiltonian_k's sign convention this session --
    running this script and comparing the printed energies against 3*t/-t/0 IS that check.

    Scope note: kite.lattice.Lattice only stores a scalar onsite energy / hopping amplitude
    per sublattice pair (see hamiltonian_k's own docstring) -- one sublattice is one orbital
    here. Fine for graphene (A, B); a genuinely multi-orbital single sublattice (as in some
    custom-vertex constructions) is out of scope for this tool.
    Last updated: 22/07/2026
"""

import numpy as np
import matplotlib.pyplot as plt

from dos_graphene import graphene_lattice
from kite import visualize as viz
import kite_style

kite_style.apply()


def main():
    t = 2.8  # eV, same convention as dos_graphene.py
    lattice = graphene_lattice(t=t)

    Gamma = np.array([0.0, 0.0])
    bz = lattice.brillouin_zone()  # 6 hexagon corners, CCW-ordered
    K = bz[0]
    M = (bz[0] + bz[1]) / 2.0  # midpoint of an edge adjacent to K

    path = viz.make_path(Gamma, K, M, Gamma, step=0.05,
                          point_labels=[r"$\Gamma$", "K", "M", r"$\Gamma$"])

    # Sanity check printed directly: these three values are what confirmed hamiltonian_k's
    # sign convention this session (see this file's own docstring).
    for label, k in (("Gamma", Gamma), ("K", K), ("M", M)):
        e = np.linalg.eigvalsh(viz.hamiltonian_k(lattice, k))
        print(f"E({label}) = {e} eV")
    print(f"Expected: E(Gamma) = +/-{3*t}, E(K) = 0, E(M) = +/-{t}")

    fig, ax = plt.subplots(figsize=(6, 5))
    viz.plot_bands(lattice, path, ax=ax, ylabel="Energy (eV)")
    ax.set_title("Graphene tight-binding bands (from the lattice definition, no KITEx)")
    fig.tight_layout()
    out_path = "plots/graphene_bands_preview.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")
    return out_path


if __name__ == "__main__":
    main()
