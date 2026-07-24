""" Spin precession in a uniform transverse field: the simplest time-evolution demo

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    The simplest possible real-time demonstration of KITE's wave-packet
    propagator: a spin-1/2 in a uniform transverse field B (a textbook Larmor
    precession problem). Two orbitals ("up","down") sit at every site of an
    ordinary square lattice with plain nearest-neighbor hopping -t (identical
    for both orbitals -- the kinetic part carries no spin information at
    all). A single on-site term couples them,

        H_onsite = B * sigma_x = [[0, B], [B, 0]],

    independent of k. Because the kinetic part is common to both orbitals and
    the field term doesn't depend on k, the Bloch Hamiltonian factorizes
    exactly as H(k) = epsilon(k)*I + B*sigma_x at every k -- so a wave packet
    seeded in the Sz=+1/2 state ("up") precesses as

        <Sz(t)> = (1/2) cos(2*B*t),   <Sy(t)> = -(1/2) sin(2*B*t),   <Sx(t)> = 0,

    EXACTLY, for any wave-packet width or shape (no dephasing, unlike a
    momentum-dependent coupling -- see the orbital-angular-momentum example
    for a case where the coupling itself depends on k).

    Operators: Sx, Sy, Sz are registered as on-site 2x2 matrices via
    add_orbital_coupling(), the same mechanism used throughout KITE's custom
    operators (custom_one, custom_two, and this generalized
    gaussian_wave_packet()) -- see docs/documentation/examples/time_evolution.md.

    Toy model, reviewed by cmt-physicist
    -------------------------------------
    B is a bare on-site coupling constant with the mathematical structure of
    a transverse Zeeman term -- there is no vector potential, no Peierls
    substitution, no g-factor, and no coupling to a real electromagnetic
    field anywhere in this construction. The hopping t is likewise pure
    numerical scaffolding (KITE's real-space method needs an extended
    periodic lattice, so it can't simulate a literal isolated two-level
    system) and carries no physics of its own here: identical for both
    orbitals, it drops out of the spin dynamics entirely. This model exists
    only to demonstrate the propagator and operator-tracking mechanism on
    the simplest Hamiltonian that shows real, exactly solvable precession.

    Last updated: 24/07/2026
"""

__all__ = ["simple_square_lattice_with_field", "main"]

import kite
import numpy as np
from kite import lattice as latt


def simple_square_lattice_with_field(t=1.0, B=0.1, a=1.0):
    """Square lattice, one site with two orbitals (up, down), a uniform
    transverse on-site field B coupling them, and plain hopping -t (same for
    both orbitals -- carries no spin information)."""
    lat = latt.Lattice(a1=[a, 0], a2=[0, a])
    pos = [0, 0]

    lat.add_sublattices(
        ('up', pos, 0.0),
        ('down', pos, 0.0),
    )
    lat.add_hoppings(
        ([1, 0], 'up', 'up', -t),
        ([0, 1], 'up', 'up', -t),
        ([1, 0], 'down', 'down', -t),
        ([0, 1], 'down', 'down', -t),
        ([0, 0], 'up', 'down', B),   # on-site transverse field, B*sigma_x
    )
    return lat


def register_spin_operators(calculation):
    """Register Sx,Sy,Sz (l0,l1,l2) as on-site 2x2 matrices."""
    calculation.add_orbital_index('up', 0)
    calculation.add_orbital_index('down', 1)

    calculation.add_orbital_coupling('down', 'up', 0.5, 'l0')
    calculation.add_orbital_coupling('up', 'down', 0.5, 'l0')

    calculation.add_orbital_coupling('down', 'up', -0.5j, 'l1')
    calculation.add_orbital_coupling('up', 'down', 0.5j, 'l1')

    calculation.add_orbital_coupling('up', 'up', 0.5, 'l2')
    calculation.add_orbital_coupling('down', 'down', -0.5, 'l2')

    return ['l0', 'l1', 'l2']


def main(t=1.0, B=0.1, num_points=400, num_moments=256, timestep=0.5,
         width=4.0, lx=128,
         output_file="spin_precession_simple-output.h5"):
    lattice = simple_square_lattice_with_field(t=t, B=B)

    nx = ny = 1
    ly = lx

    configuration = kite.Configuration(
        divisions=[nx, ny], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-5, 5],
    )
    calculation = kite.Calculation(configuration)

    operators = register_spin_operators(calculation)

    # Single k=0 wave packet, seeded purely "up" (Sz=+1/2 eigenstate).
    k_vector = np.array([[0.0, 0.0]])
    spinor = np.array([[1.0 + 0j, 0.0 + 0j]])
    mean_value = [int(lx / 2), int(ly / 2)]

    calculation.gaussian_wave_packet(
        num_points=num_points, num_moments=num_moments, num_disorder=1,
        k_vector=k_vector, spinor=spinor, width=width, timestep=timestep,
        mean_value=mean_value, operators=operators)

    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    main()
