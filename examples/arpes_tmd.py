""" Monolayer TMD, 3-band model (ARPES / spectral function)

    ##########################################################################
    #                         Copyright 2020/2022, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV, Length in nm
    Lattice: Triangular (one metal site per unit cell, 3 d-orbitals: dz2, dxy, dx2-y2)
    Configuration: Periodic boundary conditions, double precision,
                   manual scaling, size of the system 128x128, with domain decomposition (nx=ny=2)
    Calculation type: one-particle spectral function of relevance to ARPES
    Reference for the tight-binding parameters: Liu et al., Phys. Rev. B 88, 085433 (2013);
    erratum Phys. Rev. B 89, 039901 (2014).
    Last updated: 17/07/2026
"""

import kite
from kite import repository
from kite import visualize
import sys
import numpy as np


def tmd_monolayer(name="MoS2"):
    """Build a monolayer group-6 TMD lattice for the nearest-neighbor 3-band model.

    Thin wrapper around `kite.repository.group6_tmd.monolayer_3band` (`src/kite/repository.py`),
    a verified, exact port of the same nearest-neighbor 3-band tight-binding model (each of the 3
    d-orbitals -- dz2, dxy, dx2-y2 -- as its own named sublattice at the same position, connected
    by explicit scalar hoppings between every orbital pair, rather than bundling them into one
    sublattice with a dense hopping matrix). Also returns the lattice constant `a` (in nm, from
    the same parameter table), needed by `main()`'s K/Kprime corner formulas below.
    """
    lat = repository.group6_tmd.monolayer_3band(name)
    a = repository.group6_tmd._default_3band_params[name][0]
    return lat, a


def main():
    # Simulation parameters
    moments = 512
    nx = ny = 2
    lx = ly = 128

    lattice, a = tmd_monolayer("MoS2")

    # Path in k-space: K -> Gamma -> K', passing through both inequivalent valleys
    # (see the spectral-function documentation for why K and K' are inequivalent
    # here, and for the general, orientation-independent way to find them via
    # lattice.brillouin_zone()). These corner formulas are specific to the
    # a1=[a,0], a2=[a/2, a*sqrt(3)/2] convention used above.
    Gamma = np.array([0, 0])
    K = np.array([4 * np.pi / (3 * a), 0])
    Kprime = np.array([2 * np.pi / (3 * a), 2 * np.pi / (a * np.sqrt(3))])
    # `a` (and therefore K/Gamma/Kprime) is in nm here, not angstrom -- unlike this
    # script's original, locally-defined parameter table, `kite.repository.group6_tmd`
    # uses nm throughout (matching pybinding's / the rest of kite's own convention; see
    # `tmd_monolayer` above). `step` below must be expressed in the same reciprocal
    # units as the k-points (1/nm), so it is scaled up by the same factor of 10 as `a`
    # was scaled down, to keep the k-path's actual point density the same as originally
    # intended (dk = 0.1 angstrom^-1 == 1.0 nm^-1) rather than silently sampling ~10x
    # more densely (and ~10x more expensively) than before.
    dk = 1.0
    k_path = visualize.make_path(K, Gamma, Kprime, step=dk)[0]

    # One weight per orbital (3 named sublattices here, each single-orbital).
    weights = [1, 1, 1]

    # Manual spectrum range, comfortably covering the true bandwidth of this
    # model (recommended for multi-orbital models generally, rather than
    # relying on automatic scaling).
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=["periodic", "periodic"],
        is_complex=True,
        precision=1,
        spectrum_range=[-2, 4])

    calculation = kite.Calculation(configuration)
    calculation.arpes(
        k_vector=k_path,
        num_moments=moments,
        weight=weights,
        num_disorder=1)

    output_file = "arpes_tmd-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    hdf5_file = main()  # generate the Configuration file

    if len(sys.argv) > 1 and sys.argv[1] == "complete":
        sys.path.append("pybinding")  # examples/pybinding/run_all_examples.py, assumes cwd == examples/
        import run_all_examples as ra
        import process_arpes as pa
        ra.run_calculation(hdf5_file)
        ra.run_tools(hdf5_file, options="--ARPES -K jackson -S")
        pa.process_arpes("arpes.dat")
