""" Optical Conductivity of Hexagonal Boron Nitride

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Disorder: Vacancies
    Configuration: Periodic boundary conditions, double precision,
                    given rescaling, size of the system flexible, with domain decomposition (nx=ny=2)
    Calculation type: Second-order (nonlinear) optical conductivity sigma_xxy -- NOT the ordinary linear
                       sigma_ab(omega) computed by conductivity_optical elsewhere in this folder.
                       conductivity_optical_nonlinear evaluates a third-rank response tensor, which is
                       identically zero in any inversion-symmetric lattice (plain graphene, say); this
                       example uses gapped honeycomb (hexagonal boron nitride: +-gap/2 onsite potential on
                       the A/B sublattices) specifically because that staggered potential is what breaks
                       inversion symmetry and makes sigma_xxy nonzero in the first place.
                       `special=1` selects a hard-coded fast path in Src/Simulation/SimulationCondOpt2.cpp's
                       CondOpt2(): the general case computes 4 Chebyshev-moment tensors (Gamma0/1/2/3), but
                       the special=1 branch used here only computes Gamma1 and Gamma2, skipping Gamma0 and
                       Gamma3 -- per the C++ source comment, this is a materials-specific simplification for
                       HBN ("nonlinear but only has simple objects that need calculating"), not a generic
                       numerical shortcut; do not reuse special=1 for a different lattice without checking
                       whether the same simplification actually applies there.
    Last updated: 08/05/2025
"""
from kite import lattice as latt
import numpy as np
import kite


def hbn():
    t = -1.0
    a_cc = 1                              # [nm] carbon-carbon distance
    a = a_cc*np.sqrt(3)                   # [nm] unit cell length
    gap = 0.2                             # [t] gap

    lat = latt.Lattice(
        a1=[a/2,  a/2 * np.sqrt(3)],
        a2=[a/2, -a/2 * np.sqrt(3)])

    lat.add_sublattices(
    ('A', [0,    0], -gap/2.0),
    ('B', [0, a_cc],  gap/2.0)
    )

    lat.add_hoppings(
    ([ 0, 0], 'A', 'B', t),
    ([ 0, 1], 'A', 'B', t),
    ([-1, 0], 'A', 'B', t)
    )


    return lat


def main():
    lattice = hbn()

    N = 256         # number of polynomials
    C = 0.02        # concentration of vacancies
    lx = ly = 128   # system dimensions
    nx, ny = 2, 2   # number of threads in each direction
    LIM = 3.2       # overestimated bounds of the spectrum
    output_file = "hbn-output.h5"

    disorder = kite.Disorder(lattice)
    struct_A = kite.StructuralDisorder(lattice, concentration=C)
    struct_B = kite.StructuralDisorder(lattice, concentration=C)
    struct_A.add_vacancy('A')
    struct_B.add_vacancy('B')
    disorder_structural = [struct_A, struct_B]


    configuration = kite.Configuration(
            divisions=[nx, ny],
            length=[lx, ly],
            boundaries=["periodic", "periodic"],
            is_complex=True,
            precision=1,
            spectrum_range=[-LIM,LIM])

    calculation = kite.Calculation(configuration)
    calculation.conductivity_optical_nonlinear(
            "xxy",
            num_points=512,
            num_moments=N,
            num_random=1,
            num_disorder=1,
            special=1
            )

    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder_structural=disorder_structural)
    return output_file

if __name__ == "__main__":
    main()
