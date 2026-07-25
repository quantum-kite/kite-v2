""" Cross-validation of the operator-weighted local spectral density.

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    KITE's LDOS methods -- calculation.ldos() (exact/deterministic Chebyshev
    moments) and calculation.ldos_map() (stochastic/Markov full-lattice map)
    -- both accept an operators= argument that generalizes the plain local
    density of states LDOS(r,E) = -1/pi*Tr[Im G(r,r,E)] to an operator-
    weighted local spectral density Tr[O_r*Im G(r,r,E)], for any on-site
    (single-site, possibly orbital-mixing) Hermitian operator O_r registered
    via add_orbital_index()/add_orbital_coupling(). This is the real-space,
    energy-resolved analogue of a spin-resolved ARPES measurement: the local
    expectation value of an operator AT A FIXED ENERGY, not integrated over
    the Fermi function.

    This script runs BOTH methods on the same small system and the same
    operator, so the results can be cross-checked directly against each
    other -- the exact method has no statistical error at all (a Chebyshev
    moment calculation), while the stochastic method converges to it only in
    the vectors_ -> infinity limit, so agreement is checked within the
    stochastic method's own reported error bars, not exactly.

    System: the same 2D d-wave altermagnet used in altermagnet_arpes.py
    (examples/altermagnet_arpes.py) -- two magnetic sublattices A, B with a
    checkerboard arrangement, spin-doubled into 4 sublattices Aup, Bup, Adn,
    Bdn. It has zero net magnetic moment (Tr H(k) = 0 at every k) but is
    spin-split at any single k -- exactly the kind of system where a plain
    LDOS shows nothing spin-specific but an Sz-weighted local spectral
    density reveals the up/down asymmetry directly, site by site.

    Two operators are registered:
      - l0 = Sz: diag(+0.5, +0.5, -0.5, -0.5) on (Aup, Bup, Adn, Bdn) --
        the physically meaningful quantity.
      - l1 = a single-orbital projector onto Aup: diag(1, 0, 0, 0) -- a
        regression check. Since Tr[P*Im G] for a single-orbital projector P
        is mathematically identical to the plain per-orbital LDOS at that
        orbital, both operator-weighted outputs must equal (exact method)
        or statistically agree with (stochastic method) the corresponding
        slice of the plain (operators=None) output. This is the concrete
        sanity check that the operator machinery reduces correctly to the
        pre-existing diagonal LDOS in the one case where they must agree.

    Why a magnetic impurity is needed for a nonzero Sz(r)
    ------------------------------------------------------
    A CLEAN altermagnet has Sz(r,E) = 0 EXACTLY at every site and energy: the
    model's defining symmetry is (90-degree real-space rotation) x (spin
    swap up<->down), and this combined symmetry alone already forces the
    up/down contributions to cancel locally, not just after summing over the
    Brillouin zone. A generic non-magnetic defect (e.g. a vacancy) that
    itself respects that combined symmetry would leave Sz(r) = 0 even right
    at the defect -- it doesn't couple to spin, so it can't break the part
    of the symmetry that enforces the cancellation.

    To see a genuine real-space Sz(r) texture, this script instead adds a
    single MAGNETIC impurity: a strong onsite potential shift applied to
    ONLY the 'Aup' orbital (not 'Adn') at one site, via
    kite.StructuralDisorder. This directly breaks the spin-swap part of the
    symmetry at that site, so Sz polarizes locally and (Friedel-oscillation
    style) decays away from it -- exactly the kind of spatial structure a
    spin-resolved local probe is meant to reveal. calculation.ldos() is
    requested at four positions at increasing distance from the impurity
    along x, and calculation.ldos_map() maps Sz(r) over the whole lattice at
    the same energy, so the decay is visible directly in the map.

    Last updated: 24/07/2026
"""

import sys
import os

import kite
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from altermagnet_arpes import altermagnet

__all__ = ["main"]


def register_operators(calculation):
    """Register Sz (l0) and a single-orbital Aup-projector (l1).

    Sublattice order fixes the orbital index: Aup=0, Bup=1, Adn=2, Bdn=3
    (see altermagnet()'s docstring in altermagnet_arpes.py).
    """
    for label, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(label, idx)

    calculation.add_orbital_coupling('Aup', 'Aup', 0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bup', 0.5, 'l0')
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l0')

    calculation.add_orbital_coupling('Aup', 'Aup', 1.0, 'l1')

    return ['l0', 'l1']


def main(t=1.0, J=0.8, delta=0.4, lx=32, ly=32, limp=-3.0,
         energies=(-2.0, 0.0, 2.0), num_moments=512, num_disorder=1,
         sigma=0.3, vectors=4000,
         output_file="ldos_spin_operator_validation-output.h5"):
    lattice = altermagnet(t, J, delta)

    # Single magnetic impurity at the lattice center: a strong onsite
    # potential shift on 'Aup' only (not 'Adn') -- see the module docstring
    # for why this, and not a plain vacancy, is needed to get a nonzero
    # real-space Sz(r) texture out of this model.
    impurity_pos = [lx // 2, ly // 2]
    struc_disorder = kite.StructuralDisorder(lattice, position=[impurity_pos])
    struc_disorder.add_structural_disorder((impurity_pos, 'Aup', limp))

    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-6, 6],
    )
    calculation = kite.Calculation(configuration)

    operators = register_operators(calculation)

    # Four positions at increasing distance from the impurity along x (all
    # on the impurity's own sublattice, Aup, so the sequence isolates the
    # spatial decay rather than mixing in the trivial Aup/Adn offset).
    ix, iy = impurity_pos
    distances = [0, 1, 4, 8]
    positions = [[(ix + d) % lx, iy] for d in distances]

    calculation.ldos(
        energy=list(energies), num_moments=num_moments,
        position=positions, sublattice=['Aup'],
        num_disorder=num_disorder, operators=operators)

    calculation.ldos_map(
        energy_=energies[len(energies) // 2], sigma_=sigma, vectors_=vectors,
        operators=operators)

    kite.config_system(lattice, configuration, calculation,
                        disorder_structural=struc_disorder, filename=output_file)
    return output_file


if __name__ == "__main__":
    fname = main()
    print(f"Wrote {fname}")
    print("Run:   ../build/KITEx", fname)
    print("Then:  ../build/KITE-tools", fname, "--LDOS -K jackson")
    print("Then:  python ldos_spin_operator_validation_process.py", fname)
