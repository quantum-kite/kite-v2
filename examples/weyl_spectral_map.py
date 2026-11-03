""" Weyl semimetal: momentum-space spectral-function map at fixed energy

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Two-sublattice (A, B) cubic lattice with a mass term (+/- m*t onsite) and
    complex nearest-neighbour hoppings along x and y (odd powers of i on the
    A-B / B-A bonds) plus a real hopping along z. This combination realizes a
    pair of Weyl nodes in the 3D band structure -- point-like band touchings
    with linear (Dirac-like) dispersion in all three momentum directions.

    calculation.spectral_map computes the momentum-resolved spectral function
    A(k, E) on the reciprocal-space grid dual to the sample's real-space
    supercell, at a fixed energy E=0.7 (away from the Weyl node itself, so it
    probes the linearly-dispersing cone rather than the touching point).
    boundaries=["open", "random", "random"] leaves x as an open (real-space,
    slab-like) direction and resolves ky, kz as (twist-averaged) momenta --
    this is the paper authors' own intended usage for spectral_map, so a
    single spectral_map output is naturally a (ky, kz) slice at each x-layer;
    see main()'s comments and the analysis snippet for how this is read out.

    Method / reference
    -------------------
    spectral_map is a linear-scaling, stochastic (random-vector) estimator of
    the momentum-resolved spectral function at every k-point simultaneously.
    `vectors_` controls the number of random-vector realizations averaged
    over; per-k-point sampling error from a finite `vectors_` is rigorously
    bounded via a Markov-inequality argument, per:

        H. P. Veiga, D. R. Pinheiro, J. P. Santos Pires,
        J. M. Viana Parente Lopes, "Markov Inequality as a Tool for
        Linear-Scaling Estimation of Local Observables," arXiv:2510.21688,
        Phys. Rev. Research (doi 10.1103/qb1w-44r1).

    Units: energy in units of hopping |t| = 1.
    Sizing note: shrunk from a research-scale 128^3 lattice down to a
    laptop-runnable 48^3 (see main()) to keep wall-clock in the
    seconds-to-minutes range; the method is linear-scaling and stochastic, so
    it tolerates this reduction well (just a somewhat noisier map).
    Last updated: 22/07/2026
"""

__all__ = ["weyl_semimetal", "main"]

import sys

import numpy as np

import kite
from kite import lattice as latt


def weyl_semimetal(t=1.0, m=2.0):
    """Return the two-sublattice cubic lattice hosting a pair of Weyl nodes."""
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    lat = latt.Lattice(a1=a1, a2=a2, a3=a3)
    lat.add_sublattices(("A", [0, 0, 0], m * t), ("B", [0, 0, 0], -m * t))
    lat.add_hoppings(
        ([1, 0, 0], "A", "A", -0.5 * t),
        ([1, 0, 0], "B", "B", 0.5 * t),
        ([1, 0, 0], "A", "B", -0.5j * t),
        ([1, 0, 0], "B", "A", -0.5j * t),
        ([0, 1, 0], "A", "A", -0.5 * t),
        ([0, 1, 0], "B", "B", 0.5 * t),
        ([0, 1, 0], "A", "B", -0.5 * t),
        ([0, 1, 0], "B", "A", 0.5 * t),
        ([0, 0, 1], "A", "A", -0.5 * t),
        ([0, 0, 1], "B", "B", 0.5 * t),
    )
    return lat


def main(t=1.0, m=2.0, output_file="weyl_spectral_map-output.h5"):
    """Prepare the input file for KITEx (momentum-space spectral map at E=0.7)."""
    lattice = weyl_semimetal(t, m)

    # domain-decomposition parts [nx, ny, nz] (parallelization only, not
    # physics). nx = ny = nz = 2 (the original script had a typo here:
    # `nx = ny = nx = 2` assigned nx twice and left nz undefined).
    nx = ny = nz = 2

    # number of unit cells in each direction. Shrunk from the original
    # 128x128x128 (research scale, ~4.2M sites x 128 vectors) down to
    # 48x48x48, which still resolves a smooth (ky, kz) spectral map while
    # running in well under a minute on a laptop.
    lx = ly = lz = 48

    mode = "random"
    configuration = kite.Configuration(
        divisions=[nx, ny, nz],
        length=[lx, ly, lz],
        boundaries=["open", mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-6.0, 6.0],
    )
    calculation = kite.Calculation(configuration)
    # vectors_ shrunk from 128 to 32 -- still enough random-vector samples to
    # resolve the spectral map, just a bit noisier (expected and acceptable
    # per the Markov-inequality error bound, see docstring).
    calculation.spectral_map(energy_=0.7, sigma_=1e-2, vectors_=32)

    # kept in examples/ (no "Data/" prefix, unlike the original) so the
    # script runs standalone without needing to pre-create a subdirectory.
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the output, run:
    # ../build/KITEx weyl_spectral_map-output.h5
    # There is no KITE-tools post-processing step for spectral_map -- read
    # the /Calculation/spectral_map/Map dataset directly from the output
    # HDF5 with h5py.
    return output_file


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
