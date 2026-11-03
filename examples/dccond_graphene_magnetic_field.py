""" Graphene in a perpendicular magnetic field: quantum Hall conductivity plateaus

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping |t| = 1, length in units of lattice parameter |a| = 1
           (dimensionless -- see note below on why real nm/eV units cannot be used here)
    Lattice: Honeycomb (same lattice/hoppings as dos_graphene.py, just dimensionless a=t=1)
    Configuration: Periodic boundary conditions, size of the system 128x128 (nx=ny=2),
                   real Peierls-substituted magnetic field via kite.Modification
    Disorder: Small uniform on-site disorder on both sublattices (helps KPM convergence;
              kept well below the expected Landau-level gaps so plateaus stay flat)
    Calculation type: DOS and DC (xy) Hall conductivity
    Last updated: 22/07/2026

    Why dimensionless units, not real nm/eV (like dos_graphene.py)
    ---------------------------------------------------------------
    kite.Modification(magnetic_field=...) quantizes the requested field to the nearest
    integer multiple of a minimum field set by the sample's real physical area
    (Src/kite/__init__.py's config_system: B_min = Phi_0 / (ly * unit_cell_area)). With
    graphene's REAL lattice constant (~0.246 nm) at lx=ly=128, B_min alone is already
    several hundred Tesla -- physically absurd, and no laptop-scale sample can do better.
    Using dimensionless a=1 instead (same convention as dos_dccond_square_lattice.py,
    which uses magnetic_field=40 on exactly this basis) makes B_min an O(10-100) number
    in the code's internal units, so a modest-looking "magnetic_field" argument can
    correspond to a physically sensible number of flux quanta. The physics (dispersion
    shape, plateau existence) is governed by dimensionless ratios (flux quanta count,
    hopping-to-disorder ratio), so the choice of length/energy unit is immaterial here.

    Parameter choices (fixed, not auto-tuned -- see docstring of main())
    ---------------------------------------------------------------------
    lx=ly=128 (fixed by explicit request, not chosen for convenience). With this size,
    B_min = Phi_0 / (128 * unit_cell_area) ~= 37.3 (internal units); magnetic_field is
    chosen close to 10 * B_min, i.e. ~10 flux quanta threading the sample -- enough
    states per Landau level for reasonably converged statistics while keeping the level
    spacing coarse enough that only a handful of levels sit inside a modest energy window
    (the goal here is 5-7 resolved levels total, not a large Landau fan).

    num_moments=512 is deliberately modest -- this example is not trying to resolve fine
    structure between closely-spaced high-index Landau levels, just the first few on each
    side of the Dirac point. The harder convergence problem for a stochastic-trace
    quantity like conductivity_dc is NOT moment count but random-vector sampling noise:
    following the precedent in dos_dccond_haldane.py (whose conductivity_dc call uses
    num_random=10 while its DOS call only needs num_random=1), this example uses
    num_random=15 for conductivity_dc specifically.
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
