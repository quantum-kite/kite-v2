""" Kane-Mele model: spin Hall effect via KITE custom vertex operators

    ##########################################################################
    #                         Copyright 2020/2025, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    The Kane-Mele model is the honeycomb-lattice Z2 topological insulator: two
    time-reversed copies of the Haldane model, one per spin. Intrinsic
    spin-orbit coupling opens a bulk gap ~ 6*sqrt(3)*t2 and drives a quantized
    spin Hall response. Here we compute the spin Hall conductivity
    sigma^{s_z}_{xy} directly, using KITE's rank-two custom-vertex (Kubo-Bastin)
    machinery, and look for the expected FLAT PLATEAU of the response across the
    bulk gap.

    The two operators fed to custom_two are:
      A = (1/2){v_x, s_z}   (the spin current in x, symmetrized anticommutator)
      B = v_y               (the charge velocity in y)
    KITE then builds the double-Chebyshev moment matrix Gamma_{mn} =
    Tr[T_m(H) B T_n(H) A]; the Kubo-Bastin energy integral is done in Python
    (see kane_mele_spin_hall_process.py, adapted from examples/cond_sum.py, the
    diag=0 antisymmetric/Hall branch).

    Sign convention
    ---------------
    The s_z operator ("l0") is diagonal with eigenvalues +1/2 on the spin-up
    sublattices (Aup, Bup) and -1/2 on the spin-down sublattices (Adn, Bdn),
    i.e. s_z = (1/2) sigma_z (hbar = 1). The OVERALL SIGN of the computed spin
    Hall response is fixed by this choice: flipping the sign of every l0 matrix
    element flips the sign of the plateau. The Haldane phase for spin up is
    +i (t2 * (-1j) on A-A bonds, following the reference Haldane example),
    and spin down carries the opposite phase -- this opposite SOC sign between
    spins is what makes this Kane-Mele (a Z2 insulator with a spin Hall
    plateau) rather than two stacked Haldane copies (a charge Chern insulator).

    Operator-string convention
    ---------------------------
    Operator strings apply RIGHT-TO-LEFT: "vx.l0" = v_x . l0, with l0 acting
    first. The symmetrized spin current (1/2){v_x, s_z} therefore needs BOTH
    orderings, [[0.5,"vx.l0"],[0.5,"l0.vx"]]. num_moments must be EVEN.

    Every one of the 4 sublattices is registered with add_orbital_index (the
    orbital operator is a full N_orb x N_orb matrix; registering a subset would
    silently truncate it). is_complex=True is MANDATORY: the imaginary SOC
    hoppings and the (implicitly real here, but the machinery is complex) l0
    matrix would otherwise be silently dropped.

    Units: energy in units of hopping |t| = 1, length in nm.
    Last updated: 17/07/2026
"""

__all__ = ["kane_mele", "main"]

import kite
import numpy as np
from kite import custom
from kite import lattice as latt


def kane_mele(t=1.0, t2=0.15):
    """Return the spin-full Kane-Mele honeycomb lattice (4 sublattices).

    Sublattice add order fixes alias_id: Aup=0, Bup=1, Adn=2, Bdn=3.
    The up-spin block reproduces the Haldane model of
    examples/dos_dccond_haldane.py exactly; the down-spin block is its complex
    conjugate (opposite Haldane phase).
    """
    a = 0.24595   # [nm] unit-cell length
    a_cc = 0.142  # [nm] carbon-carbon distance

    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, 0.5 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)

    # Order matters: alias_id 0,1,2,3
    lat.add_sublattices(
        ('Aup', [0, -a_cc / 2], 0.0),
        ('Bup', [0,  a_cc / 2], 0.0),
        ('Adn', [0, -a_cc / 2], 0.0),
        ('Bdn', [0,  a_cc / 2], 0.0),
    )

    lat.add_hoppings(
        # --- Nearest-neighbour, spin-diagonal (charge hopping), both spins ---
        ([0, 0], 'Aup', 'Bup', -t),
        ([1, -1], 'Aup', 'Bup', -t),
        ([0, -1], 'Aup', 'Bup', -t),
        ([0, 0], 'Adn', 'Bdn', -t),
        ([1, -1], 'Adn', 'Bdn', -t),
        ([0, -1], 'Adn', 'Bdn', -t),

        # --- Next-nearest-neighbour intrinsic SOC, spin up (Haldane phase +i) ---
        ([1, 0], 'Aup', 'Aup', -t2 * 1j),
        ([0, -1], 'Aup', 'Aup', -t2 * 1j),
        ([-1, 1], 'Aup', 'Aup', -t2 * 1j),
        ([1, 0], 'Bup', 'Bup', -t2 * -1j),
        ([0, -1], 'Bup', 'Bup', -t2 * -1j),
        ([-1, 1], 'Bup', 'Bup', -t2 * -1j),

        # --- NNN SOC, spin down: OPPOSITE imaginary sign (complex conjugate) ---
        ([1, 0], 'Adn', 'Adn', -t2 * -1j),
        ([0, -1], 'Adn', 'Adn', -t2 * -1j),
        ([-1, 1], 'Adn', 'Adn', -t2 * -1j),
        ([1, 0], 'Bdn', 'Bdn', -t2 * 1j),
        ([0, -1], 'Bdn', 'Bdn', -t2 * 1j),
        ([-1, 1], 'Bdn', 'Bdn', -t2 * 1j),
    )
    return lat


def main(t=1.0, t2=0.15, moments=256, output_file="kane_mele_spin_hall-output.h5"):
    """Prepare the input file for KITEx (spin Hall via custom_two)."""
    lattice = kane_mele(t, t2)

    nx = ny = 2
    lx = ly = 64

    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=['periodic', 'periodic'],
        is_complex=True,
        precision=1,
        spectrum_range=[-4, 4],
    )

    calculation = kite.Calculation(configuration)

    # Register EVERY sublattice (idx = add order = alias_id).
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)

    # s_z = (1/2) sigma_z: +1/2 on spin up, -1/2 on spin down.
    calculation.add_orbital_coupling('Aup', 'Aup',  0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bup',  0.5, 'l0')
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l0')

    # A = (1/2){v_x, s_z} spin current ; B = v_y charge velocity
    A = custom.Vertex(moments, [[0.5, "vx.l0"], [0.5, "l0.vx"]])
    B = custom.Vertex(moments, [[1.0, "vy"]])

    calculation.custom_two(
        stream_=[A, B],
        num_random_=1,
        num_disorder_=1,
        num_points_=1000,
        temperature_=0.01,
    )

    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # Run:   ../build/KITEx kane_mele_spin_hall-output.h5
    # There is NO KITE-tools step for custom_two -- the Kubo-Bastin energy
    # integral is done in Python; see kane_mele_spin_hall_process.py.
    return output_file


if __name__ == "__main__":
    main()
