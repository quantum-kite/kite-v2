""" Real-time precession of orbital angular momentum and orbital quadrupoles

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Orbital quadrupoles Q_ab = {L_a, L_b}/2 describe the "shape" of an
    orbital state, one rank higher than the orbital angular momentum L
    itself. Under time evolution the two are coupled: the Heisenberg torque
    on Q_ab, T_Qab = (i/hbar)[H, Q_ab], is generally nonzero whenever the
    Hamiltonian's hopping is anisotropic between orbital characters (here,
    V_pp_sigma != V_pp_pi) -- that anisotropy is exactly what breaks the
    rotational symmetry that would otherwise decouple L and Q. This example
    tracks <L(t)> and <Q(t)> simultaneously on the same real-time wave-packet
    propagation to show that coupling directly.

    Model: same sp3 Slater-Koster square lattice (s, p_x, p_y, p_z on-site
    orbitals) used for the orbital Hall effect example in this folder, with
    spin-orbit coupling switched off -- the p_x, p_y, p_z manifold is exactly
    the l=1 angular-momentum representation needed to build L and Q.

    Operators: L_x, L_y, L_z and the six Q_ab are all built as on-site
    orbital-space matrices via add_orbital_coupling(), the same registration
    mechanism (and the same on-site-only restriction -- see
    docs/api/kite.md's warning under add_orbital_coupling()) used for spin in
    gaussian_wavepacket_only_anderson.py. Convention: (L_k)_{jl} = -i *
    epsilon_{kjl} on the (p_x, p_y, p_z) block (verified: Hermitian,
    [L_i, L_j] = i * epsilon_ijk L_k, L^2 = 2*1 -- the standard l=1 algebra).
    This differs by an overall sign from the L_z-like "l0" operator used in
    the orbital Hall example; as with spin, the overall sign of an angular
    momentum operator is a convention, not a physical ambiguity.

    Validation (reviewed by cmt-physicist)
    ---------------------------------------
    With the default initial spinor below (a real p_x,p_y combination plus a
    small p_z admixture, no p_x - i*p_y component), <Q_xy(t)> is exactly zero
    at every timestep, not just at t=0 -- a genuine selection rule, not a
    bug. The lattice's discrete C4 rotation (px->py, py->-px, pz->pz,
    x<->y bonds swapped) commutes with H even though V_pp_sigma != V_pp_pi
    breaks the full continuous rotational symmetry; the initial state only
    populates the C4 eigenvalue sectors connected to m=+1 and m=0 (never
    m=-1), and Q_xy only couples sectors differing by |Delta m|=2 -- so it
    stays unpopulated for all time. Q_xz and Q_yz (which couple |Delta m|=1,
    populated here) are generically nonzero and unequal to each other. This
    is a cheap, checkable prediction of the operator construction, not an
    accidental cancellation.

    Last updated: 23/07/2026
"""

__all__ = ["sp3_square_lattice", "main"]

import kite
import numpy as np
from kite import lattice as latt

# Slater-Koster sp3 parameters (SP3-Go from sp3-model-OHE.py; SOC switched off).
a = 1.0
E_s, E_px, E_py, E_pz = 3.2, -0.5, -0.5, -0.5
V_ss_sigma = -0.5
V_pp_sigma = 0.5
V_pp_pi = -0.2
V_sp_sigma = 0.5


def sp3_square_lattice():
    """Square lattice, one atom per cell, s/p_x/p_y/p_z orbitals (Slater-Koster)."""
    lat = latt.Lattice(a1=[a, 0], a2=[0, a])
    pos = [0, 0]

    lat.add_sublattices(
        ('s', pos, E_s),
        ('px', pos, E_px),
        ('py', pos, E_py),
        ('pz', pos, E_pz),
    )

    lat.add_hoppings(
        ([1, 0], 's', 's', V_ss_sigma),
        ([1, 0], 's', 'px', V_sp_sigma),
        ([1, 0], 'px', 'px', V_pp_sigma),
        ([1, 0], 'py', 'py', V_pp_pi),
        ([1, 0], 'pz', 'pz', V_pp_pi),
        ([-1, 0], 's', 'px', -V_sp_sigma),
    )
    lat.add_hoppings(
        ([0, 1], 's', 's', V_ss_sigma),
        ([0, 1], 's', 'py', V_sp_sigma),
        ([0, 1], 'px', 'px', V_pp_pi),
        ([0, 1], 'py', 'py', V_pp_sigma),
        ([0, 1], 'pz', 'pz', V_pp_pi),
        ([0, -1], 's', 'py', -V_sp_sigma),
    )
    return lat


def register_L_and_Q(calculation):
    """Register L_x,L_y,L_z (l0,l1,l2) and the six Q_ab (l3..l8) as on-site operators."""
    calculation.add_orbital_index('s', 0)
    calculation.add_orbital_index('px', 1)
    calculation.add_orbital_index('py', 2)
    calculation.add_orbital_index('pz', 3)

    # L_x
    calculation.add_orbital_coupling('pz', 'py', -1j, 'l0')
    calculation.add_orbital_coupling('py', 'pz', 1j, 'l0')
    # L_y
    calculation.add_orbital_coupling('px', 'pz', -1j, 'l1')
    calculation.add_orbital_coupling('pz', 'px', 1j, 'l1')
    # L_z
    calculation.add_orbital_coupling('py', 'px', -1j, 'l2')
    calculation.add_orbital_coupling('px', 'py', 1j, 'l2')

    # Q_xx = L_x^2 = diag(0, 1, 1) on (px, py, pz)
    calculation.add_orbital_coupling('py', 'py', 1.0, 'l3')
    calculation.add_orbital_coupling('pz', 'pz', 1.0, 'l3')
    # Q_yy = diag(1, 0, 1)
    calculation.add_orbital_coupling('px', 'px', 1.0, 'l4')
    calculation.add_orbital_coupling('pz', 'pz', 1.0, 'l4')
    # Q_zz = diag(1, 1, 0)
    calculation.add_orbital_coupling('px', 'px', 1.0, 'l5')
    calculation.add_orbital_coupling('py', 'py', 1.0, 'l5')
    # Q_xy = {L_x,L_y}/2
    calculation.add_orbital_coupling('py', 'px', -0.5, 'l6')
    calculation.add_orbital_coupling('px', 'py', -0.5, 'l6')
    # Q_xz = {L_x,L_z}/2
    calculation.add_orbital_coupling('pz', 'px', -0.5, 'l7')
    calculation.add_orbital_coupling('px', 'pz', -0.5, 'l7')
    # Q_yz = {L_y,L_z}/2
    calculation.add_orbital_coupling('pz', 'py', -0.5, 'l8')
    calculation.add_orbital_coupling('py', 'pz', -0.5, 'l8')

    return ['l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8']


def main(num_points=400, num_moments=512, timestep=0.05, width=2.0,
         output_file="oam_quadrupole_precession-output.h5"):
    lattice = sp3_square_lattice()

    nx = ny = 1
    lx = ly = 16

    configuration = kite.Configuration(
        divisions=[nx, ny], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-8, 8],
    )
    calculation = kite.Calculation(configuration)

    operators = register_L_and_Q(calculation)

    # Single k=0 wave packet, initial state with sizeable L_z (p_x + i*p_y-like
    # combination) so the L-Q coupling has something to precess from.
    k_vector = np.array([[0.0, 0.0]])
    spinor = np.array([[0.0, 1.0 + 0j, 1j, 0.3 + 0j]])
    spinor /= np.linalg.norm(spinor)
    mean_value = [int(lx / 2), int(ly / 2)]

    calculation.gaussian_wave_packet(
        num_points=num_points, num_moments=num_moments, num_disorder=1,
        k_vector=k_vector, spinor=spinor, width=width, timestep=timestep,
        mean_value=mean_value, operators=operators)

    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    main()
