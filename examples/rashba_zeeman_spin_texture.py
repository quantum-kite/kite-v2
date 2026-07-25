""" Real-space spin texture (Sx, Sy, Sz) around a vacancy in a Rashba+Zeeman
    honeycomb lattice.

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Reuses the spin-doubled honeycomb lattice (Aup, Bup, Adn, Bdn) and Rashba
    SOC hopping construction from examples/rashba_edelstein_graphene.py
    unchanged, adding an out-of-plane Zeeman splitting Bz (onsite +Bz on
    Aup/Bup, -Bz on Adn/Bdn) -- the standard "Rashba + Zeeman" model, the same
    combination used by Qiao et al., "Quantum Anomalous Hall Effect in
    Graphene from Rashba and Exchange Effects", Phys. Rev. B 82, 161414(R)
    (2010) [arXiv:1005.1672], to open a topological gap in honeycomb graphene
    (there studied for its Chern number/edge states; here for the resulting
    local spin texture around a defect instead).

    Why a vacancy (not a magnetic impurity) should work here, unlike the
    altermagnet case
    ------------------------------------------------------------------------
    In examples/ldos_spin_operator_validation.py, a clean altermagnet has
    Sz(r)=0 EXACTLY everywhere because of a combined (rotation)x(spin-swap)
    symmetry, and a plain (spin-independent) vacancy was shown NOT to break
    it. Rashba SOC is different: it locks spin direction to real-space
    hopping direction (the vertex operator itself, e.g.
    coefficient ~ (dy +/- i*dx), depends on the bond's spatial orientation),
    so spin is not a separate, decoupled label the way it is in the
    altermagnet model -- ANY real-space scatterer, magnetic or not, couples
    to spin through this locking. A single nonmagnetic vacancy is therefore
    expected to induce a genuine spin density texture around it (a real-space
    analogue of the well-documented spin-Hall/spin-swapping response of
    Rashba two-dimensional electron gases to nonmagnetic defects) -- this is
    the hypothesis this script tests, not an established KITE result.

    Calculation
    -----------
    calculation.ldos_map(operators=['l0','l1','l2']) maps Sx, Sy, Sz
    (registered via add_orbital_coupling, same mechanism as
    rashba_edelstein_graphene.py's l0) over the whole lattice at one energy
    inside the band. Post-process with rashba_zeeman_spin_texture_process.py.

    Units: energy in units of hopping |t|=1, length in nm.
    Last updated: 25/07/2026
"""

import sys
import os

import kite
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from rashba_edelstein_graphene import rashba_graphene

__all__ = ["rashba_zeeman_honeycomb", "register_spin_operators", "main"]


def rashba_zeeman_honeycomb(t=1.0, t2=0.0, lambda_R=0.3, Bz=0.3):
    """rashba_graphene() plus an out-of-plane Zeeman onsite splitting Bz."""
    lat = rashba_graphene(t, t2, lambda_R)
    for name, sub in lat.sublattices.items():
        sign = +1.0 if name in ('Aup', 'Bup') else -1.0
        sub.energy = np.array([[sign * Bz]])
    return lat


def register_spin_operators(calculation):
    """Register Sx, Sy, Sz SEPARATELY on the A and B sublattices (6 operators
    total: l0..l2 = Sx_A,Sy_A,Sz_A; l3..l5 = Sx_B,Sy_B,Sz_B).

    ldos_map's operator-weighted output is one value per unit cell, summed
    over whichever orbitals the requested operator has nonzero entries on
    (see docs/api/kite.md's operators= entry). A single 4-orbital Sx/Sy/Sz
    matrix spanning BOTH A and B therefore collapses the two-atom honeycomb
    basis into one number per unit cell -- exactly the (visually
    unsatisfying, symmetry-hiding) triangular-lattice-only result from the
    first version of this script. Restricting each matrix's nonzero entries
    to ONLY the A block (Aup/Adn) or ONLY the B block (Bup/Bdn) instead
    isolates each sublattice's own spin density at the SAME per-unit-cell
    resolution, so it can be plotted at that atom's own real position
    (A at the unit cell origin, B offset by a_cc) and the actual honeycomb/
    hexagonal symmetry becomes visible.
    """
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)

    # Sx_A, Sy_A, Sz_A: (1/2)*sigma_{x,y,z} acting only on the Aup/Adn block.
    calculation.add_orbital_coupling('Adn', 'Aup', 0.5, 'l0')
    calculation.add_orbital_coupling('Aup', 'Adn', 0.5, 'l0')
    calculation.add_orbital_coupling('Adn', 'Aup', -0.5j, 'l1')
    calculation.add_orbital_coupling('Aup', 'Adn', 0.5j, 'l1')
    calculation.add_orbital_coupling('Aup', 'Aup', 0.5, 'l2')
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l2')

    # Sx_B, Sy_B, Sz_B: same, on the Bup/Bdn block.
    calculation.add_orbital_coupling('Bdn', 'Bup', 0.5, 'l3')
    calculation.add_orbital_coupling('Bup', 'Bdn', 0.5, 'l3')
    calculation.add_orbital_coupling('Bdn', 'Bup', -0.5j, 'l4')
    calculation.add_orbital_coupling('Bup', 'Bdn', 0.5j, 'l4')
    calculation.add_orbital_coupling('Bup', 'Bup', 0.5, 'l5')
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l5')

    return ['l0', 'l1', 'l2', 'l3', 'l4', 'l5']


def main(t=1.0, t2=0.0, lambda_R=0.3, Bz=0.3, lx=32, ly=32,
         energy=0.5, sigma=0.3, vectors=4000,
         output_file="rashba_zeeman_spin_texture-output.h5"):
    lattice = rashba_zeeman_honeycomb(t, t2, lambda_R, Bz)

    # Single nonmagnetic vacancy at the lattice center (Aup sublattice) --
    # see the module docstring for why this, unlike the altermagnet case,
    # should induce real-space spin texture via Rashba spin-momentum locking.
    impurity_pos = [lx // 2, ly // 2]
    struc_disorder = kite.StructuralDisorder(lattice, position=[impurity_pos])
    struc_disorder.add_vacancy('Aup')

    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-6, 6],
    )
    calculation = kite.Calculation(configuration)

    operators = register_spin_operators(calculation)

    calculation.ldos_map(
        energy_=energy, sigma_=sigma, vectors_=vectors, operators=operators)

    kite.config_system(lattice, configuration, calculation,
                        disorder_structural=struc_disorder, filename=output_file)
    return output_file


if __name__ == "__main__":
    fname = main()
    print(f"Wrote {fname}")
    print("Run:   ../build/KITEx", fname)
    print("Then:  python rashba_zeeman_spin_texture_process.py", fname)
