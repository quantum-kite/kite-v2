""" Graphene Rashba-Edelstein effect: Kane-Mele + Rashba SOC via KITE custom vertex operators

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    The Rashba-Edelstein effect (REE) is the spin-density response chi_yx =
    d<s_x>/dE_y to an in-plane electric field, present only when inversion
    symmetry is broken (Rashba SOC, lambda_R != 0). Kane-Mele intrinsic SOC
    (lambda_I) alone is inversion-symmetric and gives spin Hall, NOT REE -- it
    should contribute zero to this response. Both terms live on the same
    spin-doubled honeycomb lattice as examples/kane_mele_spin_hall.py.

    Vertex: bare spin density, not spin current --
    docs/documentation/examples/custom_vertex_operators.md's own comparison
    table already lists this exact pair for spin/orbital Edelstein responses:

        A = s_x   (bare spin operator, built via add_orbital_coupling)
        B = v_y

    Rashba hopping term: adapted, not re-derived
    -----------------------------------------------
    Rather than re-deriving the honeycomb Rashba bond phases from scratch, this
    reuses the exact numerical convention already validated and shipped
    elsewhere in this codebase -- examples/paper/Section_4_B_topology's
    qahe_disorder.py and examples/paper/Section_4_E_spintronics's
    gaussian_wavepacket_only_anderson.py both use

        rashba_so = lambda_R * 2j / 3

    with per-bond coefficients that are a pure phase times rashba_so,
    determined by the bond's unit direction d_hat=(dx,dy):

        coefficient(Aup -> Bdn, bond d_hat) = rashba_so * (dy + 1j*dx)
        coefficient(Adn -> Bup, bond d_hat) = rashba_so * (dy - 1j*dx)

    (cross-checked against BOTH reference files' explicit numbers for two
    different unit-cell orientations -- the formula reproduces every one of
    their hard-coded delta1/delta2/delta3 coefficients exactly). Those two
    files use different primitive vectors than kane_mele_spin_hall.py's
    (a1=[a,0], a2=[a/2,a*sqrt(3)/2], Aup=[0,-a_cc/2], Bup=[0,a_cc/2]), so the
    formula above -- not their literal numbers -- is what's reused; it is
    applied to kane_mele_spin_hall.py's own bond directions, computed from
    that lattice's own already-verified geometry: the three NN bonds
    (relative_index [0,0], [1,-1], [0,-1]) point at 90, -30, -150 degrees,
    i.e. d_hat in {(0,1), (sqrt(3)/2,-1/2), (-sqrt(3)/2,-1/2)}.

    Spin-flip label gotcha
    -----------------------
    Only ONE custom operator (s_x) is registered here -- unlike
    kane_mele_spin_hall.py's s_z, which coexists with nothing else. The custom
    operator label's trailing digit is a direct 0-indexed lookup into KITEx's
    vector of registered operator matrices (Src/Simulation/Custom/
    SimulationRankOne.cpp's act_with_stream reads it back with std::stoi and
    indexes directly into that vector) -- registering the only operator here
    as "l1" instead of "l0" segfaults (out-of-bounds access), it does not
    merely misname it. So it MUST be "l0" here, even though the physical
    operator is s_x, not s_z.

    Validation
    ----------
    Two structural checks, cheap and independent of any absolute-normalization
    ambiguity:
      (a) lambda_I=0, lambda_R!=0: REE response should be nonzero.
      (b) lambda_R=0, lambda_I!=0: REE response should vanish (Kane-Mele alone
          is inversion-symmetric -- no spin-charge coupling, only spin Hall).
    See rashba_edelstein_graphene_process.py.

    Unlike kane_mele_spin_hall.py (num_random_=1 is enough there -- its spin
    Hall plateau is a strong signal deep in a real bulk gap), check (b) here
    is a near-EXACT cancellation (up/down sectors are orthogonal subspaces
    when lambda_R=0, so every stochastic sample should contribute exactly
    zero to Gamma_mn) that a single random vector only satisfies on average,
    not per-sample: at num_random_=1, the residual sampling noise from check
    (b) was comparable in size to the actual signal at lambda_R=0.1 (both
    just stochastic artifacts of one random vector); num_random_=50 shrinks
    the num_random_=1 check (b) plateau roughly 20-fold while leaving the
    lambda_R=0.1 signal's plateau essentially unchanged, confirming this is
    genuine averaging convergence, not a wrong vertex/lattice.

    Units: energy in units of hopping |t|=1, length in nm.
    Last updated: 23/07/2026
"""

__all__ = ["rashba_graphene", "main"]

import kite
import numpy as np
from kite import custom
from kite import lattice as latt


def rashba_graphene(t=1.0, t2=0.15, lambda_R=0.1):
    """Return the spin-full honeycomb lattice with Kane-Mele + Rashba SOC.

    Sublattice add order fixes alias_id: Aup=0, Bup=1, Adn=2, Bdn=3 -- same
    convention as examples/kane_mele_spin_hall.py, which this reuses unchanged
    for the charge NN hopping and the Kane-Mele intrinsic SOC.
    """
    a = 0.24595   # [nm] unit-cell length
    a_cc = 0.142  # [nm] carbon-carbon distance

    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, 0.5 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)

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

        # --- Next-nearest-neighbour intrinsic (Kane-Mele) SOC, spin up ---
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

    if lambda_R != 0.0:
        rashba_so = lambda_R * 2j / 3.0
        # d_hat for each NN bond, in kane_mele_spin_hall.py's own geometry.
        d_000 = (0.0, 1.0)
        d_1m1 = (np.sqrt(3) / 2, -0.5)
        d_0m1 = (-np.sqrt(3) / 2, -0.5)

        def up_to_dn(d):
            dx, dy = d
            return rashba_so * (dy + 1j * dx)

        def dn_to_up(d):
            dx, dy = d
            return rashba_so * (dy - 1j * dx)

        lat.add_hoppings(
            ([0, 0], 'Aup', 'Bdn', up_to_dn(d_000)),
            ([1, -1], 'Aup', 'Bdn', up_to_dn(d_1m1)),
            ([0, -1], 'Aup', 'Bdn', up_to_dn(d_0m1)),

            ([0, 0], 'Adn', 'Bup', dn_to_up(d_000)),
            ([1, -1], 'Adn', 'Bup', dn_to_up(d_1m1)),
            ([0, -1], 'Adn', 'Bup', dn_to_up(d_0m1)),
        )

    return lat


def main(t=1.0, t2=0.15, lambda_R=0.1, moments=256,
         output_file="rashba_edelstein_graphene-output.h5"):
    """Prepare the input file for KITEx (Rashba-Edelstein via custom_two)."""
    lattice = rashba_graphene(t, t2, lambda_R)

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

    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)

    # s_x = (1/2)*sigma_x: off-diagonal between up/down copies of each sublattice.
    calculation.add_orbital_coupling('Adn', 'Aup', 0.5, 'l0')
    calculation.add_orbital_coupling('Aup', 'Adn', 0.5, 'l0')
    calculation.add_orbital_coupling('Bdn', 'Bup', 0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bdn', 0.5, 'l0')

    # A = bare s_x spin density ; B = v_y charge velocity
    A = custom.Vertex(moments, [[1.0, "l0"]])
    B = custom.Vertex(moments, [[1.0, "vy"]])

    calculation.custom_two(
        stream_=[A, B],
        num_random_=50,
        num_disorder_=1,
        num_points_=1000,
        temperature_=0.01,
    )

    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # Run:   ../build/KITEx rashba_edelstein_graphene-output.h5
    # There is NO KITE-tools step for custom_two -- the Kubo-Bastin energy
    # integral is done in Python; see rashba_edelstein_graphene_process.py.
    return output_file


if __name__ == "__main__":
    main()
