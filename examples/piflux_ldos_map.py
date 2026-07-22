""" Pi-flux square lattice: real-space LDOS map around a single vacancy

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Two-sublattice (A, B) square lattice with a staggered-sign hopping
    pattern (uniform -t along x, but +t on B-B and -t on A-A bonds along y).
    This staggered sign threads exactly pi flux through every plaquette,
    which folds the square-lattice band structure onto a low-energy Dirac
    cone at E=0 (same mechanism as the pi-flux/Dirac-fermion square-lattice
    models used to mimic graphene-like physics on a bipartite lattice).

    A single vacancy is placed on the A sublattice at the center of the
    sample via kite.StructuralDisorder(position=...). calculation.ldos_map
    then computes the real-space local density of states at the Dirac-point
    energy E=0, which is expected to show a localized perturbation (LDOS
    suppression exactly at the vacancy and Friedel-like oscillations decaying
    away from it) rather than a featureless uniform map.

    Method / reference
    -------------------
    ldos_map is a linear-scaling, stochastic (random-vector) estimator of the
    real-space LDOS at every site simultaneously -- it does not loop over
    sites individually. `vectors_` controls the number of random-vector
    realizations averaged over; per-site sampling error from a finite
    `vectors_` is rigorously bounded via a Markov-inequality argument, per:

        H. P. Veiga, D. R. Pinheiro, J. P. Santos Pires,
        J. M. Viana Parente Lopes, "Markov Inequality as a Tool for
        Linear-Scaling Estimation of Local Observables," arXiv:2510.21688,
        Phys. Rev. Research (doi 10.1103/qb1w-44r1).

    Units: energy in units of hopping |t| = 1.
    Sizing note: shrunk from a research-scale 2048x2048 lattice down to a
    laptop-runnable 256x256 (see main()) to keep wall-clock in the
    seconds-to-minutes range; the method is linear-scaling and stochastic, so
    it tolerates this reduction well (just a somewhat noisier map).

    Note: is_complex=True is required below even though every hopping in this
    lattice is real-valued -- ldos_map/spectral_map are only implemented for
    the complex Hamiltonian instantiation in KITE's C++ core (real hoppings
    are still handled correctly, just stored/propagated as complex numbers).
    Last updated: 22/07/2026
"""

__all__ = ["square_lattice", "main"]

import sys

import numpy as np

import kite
import myaux
from kite import lattice as latt


def square_lattice(onsite=(0, 0), t=1):
    """Return lattice specification for the pi-flux square lattice.

    Two sublattices A, B on a rectangular unit cell (a1=[2,0], a2=[0,1]).
    The alternating sign on the y-direction A-A / B-B hoppings (vs. the
    uniform sign on the x-direction A-B hoppings) is what threads pi flux
    per plaquette.
    """
    a1 = np.array([2, 0])
    a2 = np.array([0, 1])
    lat = latt.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(("A", [0, 0], onsite[0]), ("B", [1, 0], onsite[1]))
    lat.add_hoppings(
        ([0, 0], "A", "B", -t),
        ([0, 1], "B", "B", t),
        ([0, 1], "A", "A", -t),
        ([1, 0], "B", "A", -t),
    )
    return lat


def main(onsite=(0, 0), t=1, output_file="piflux_ldos_map-output.h5"):
    """Prepare the input file for KITEx (real-space LDOS map at E=0)."""
    lattice = square_lattice(onsite, t)

    # number of unit cells in each direction. Shrunk from the original
    # 2048x2048 (research scale) to 256x256, which is still comfortably
    # large enough for a single central vacancy's perturbation to decay
    # well before it reaches the boundary/wraps around, while running in
    # well under a minute on a laptop.
    lx = 256
    ly = 256

    # single vacancy at the exact center of the sample, on sublattice A.
    # myaux.mycentral(lx, ly) returns the central unit-cell index [lx//2, ly//2]
    # -- lx, ly must be defined before this call (bug in the original script).
    imp_posA = myaux.mycentral(lx, ly)
    struc_disorder_A = kite.StructuralDisorder(lattice, position=imp_posA)
    struc_disorder_A.add_vacancy("A")
    disorder_structural = [struc_disorder_A]

    # domain-decomposition parts [nx, ny] (parallelization only, not physics)
    nx = 2
    ny = 2

    mode = "random"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        # ldos_map (and spectral_map) are compiled with
        # `if constexpr (is_tt<std::complex, T>::value)` in
        # Src/Simulation/SimulationLDoS.cpp -- for a real Hamiltonian
        # (is_complex=False) that branch is compiled away entirely and the
        # calculation silently does nothing (no error, no output dataset).
        # is_complex=True is therefore MANDATORY here, even though this
        # lattice's hoppings are all real-valued.
        is_complex=True,
        precision=1,
        spectrum_range=[-3.1, 3.1],
        seed_h=1,
    )
    calculation = kite.Calculation(configuration)
    # vectors_ shrunk from 64 to 32 -- still enough random-vector samples to
    # resolve the vacancy's LDOS perturbation, just a bit noisier (expected
    # and acceptable per the Markov-inequality error bound, see docstring).
    calculation.ldos_map(energy_=0.0, sigma_=1e-2, vectors_=32)
    kite.config_system(
        lattice,
        configuration,
        calculation,
        filename=output_file,
        disorder_structural=disorder_structural,
    )

    # for generating the output, run:
    # ../build/KITEx piflux_ldos_map-output.h5
    # There is no KITE-tools post-processing step for ldos_map -- read the
    # /Calculation/ldos_map/Map dataset directly from the output HDF5 with
    # h5py; see process_piflux_ldos.py for the exact indexing convention
    # and a ready-made plotting routine.
    return output_file


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
