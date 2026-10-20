""" Graphene in a strong magnetic field: quantum Hall conductivity

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping |t| = 1, length in units of lattice parameter |a| = 1
    Lattice: Honeycomb (same lattice/hoppings as dos_graphene.py, with a = t = 1)
    Configuration: Periodic boundary conditions, size of the system 128x128 (nx=ny=2),
                   Peierls-substituted magnetic field via kite.Modification
    Disorder: Small uniform on-site disorder on both sublattices
    Calculation type: DOS and DC (xy) Hall conductivity
    Last updated: 15/10/2026

    The field is quantized to an integer multiple m of the smallest one commensurate
    with the sample, and the simulation sees the flux per unit cell, f = m / L_y. Here
    m = 10 and L_y = 128, so f = 5/64 of a flux quantum per hexagon and the magnetic
    length is 1.33 a -- the lattice (Hofstadter) regime rather than the continuum Dirac
    one. The DOS shows an isolated zeroth Landau level at E = 0 with a gap out to
    |E| = 0.83 t, and sigma_xy a single plateau at +-1 e^2/h per spin before the band
    breaks into Hofstadter subbands.
    """

__all__ = ["graphene_lattice", "main"]

import kite
import numpy as np
from kite import lattice as latt


def graphene_lattice(onsite=(0, 0), t=1.0):
    """Return the dimensionless honeycomb lattice (a=1, t=1) -- see module docstring for why."""
    a = 1.0
    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(
        ('A', [0, -0.5 / np.sqrt(3)], onsite[0]),
        ('B', [0, 0.5 / np.sqrt(3)], onsite[1]),
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', -t),
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t),
    )
    return lat


def main(t=1.0, magnetic_field=373.0, disorder_w=0.05,
         output_file="dccond_graphene_magnetic_field-output.h5"):
    """Prepare the input file for KITEx (DOS + DC Hall conductivity in a magnetic field)."""
    lattice = graphene_lattice(t=t)

    disorder = kite.Disorder(lattice)
    disorder.add_disorder('A', 'Uniform', 0.0, disorder_w)
    disorder.add_disorder('B', 'Uniform', 0.0, disorder_w)

    nx = ny = 2
    lx = ly = 128

    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=['periodic', 'periodic'],
        is_complex=True,
        precision=1,
        spectrum_range=[-4.0, 4.0],
    )
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1000,
        num_moments=512,
        num_random=1,
        num_disorder=1,
    )
    calculation.conductivity_dc(
        num_points=1000,
        num_moments=512,
        num_random=15,
        num_disorder=1,
        direction='xy',
        temperature=0.05,
    )

    modification = kite.Modification(magnetic_field=magnetic_field)
    kite.config_system(lattice, configuration, calculation, modification=modification,
                        filename=output_file, disorder=disorder)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx dccond_graphene_magnetic_field-output.h5
    # ../build/KITE-tools dccond_graphene_magnetic_field-output.h5 --DOS --CondDC -F -4 4 400

    return output_file


if __name__ == "__main__":
    main()
