""" Real-space Sz(r) map of a bearded/Klein-terminated graphene ribbon with a
    fixed (not self-consistently solved) staggered Neel/AFM mass term, clean
    vs. a low concentration of real vacancies.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    See the module docstring of afm_zigzag_bands.py for the full physical
    picture (bulk gap at both valleys, near-flat edge bands close to
    E=+/-Delta over most of the zone -- exact at k=0, only visibly dispersive
    near where they merge with the bulk continuum -- each essentially one
    sublattice/spin at one edge, yet globally spin-up DOS == spin-down DOS by
    symmetry between the two edges). This script computes the corresponding
    local, real-space-resolved probe: calculation.ldos_map(operators=['l0','l1'])
    maps Sz_A(r), Sz_B(r) separately (see register_sz_operators()).

    Energy: the default E=0.05t probes broadened in-gap spectral weight,
    below the edge-band energy E=Delta=0.3t rather than on it. At E=Delta,
    the Gaussian window also overlaps bulk states near the band edge, making
    a single-disorder-realization comparison less visually selective -- see
    the in-depth write-up (docs/documentation/examples/afm_zigzag_ribbon.md)
    for the full derivation. For:
      - a clean ribbon (edge polarization visible, cancels in the global
        DOS but not site-by-site);
      - the same ribbon with a low (~5%) concentration of real vacancies
        (both spin channels removed at the same site -- see
        register_vacancies()) to show vacancy-induced local spin imbalance,
        especially near the edges.

    Note on terminology: the output is an energy-resolved, operator-weighted
    local spectral density (Tr[Sz * Im G(r,r,E)] at a fixed E), not an
    equilibrium spin density -- the latter would require occupation weighting
    and integration over energy, which this calculation does not do.

    Boundary conditions: periodic along a1 (armchair direction is not used
    here), open along a2 -- this produces two bearded/Klein-terminated edges
    (row y=0 and row y=Ly-1; see afm_zigzag_bands.py's module docstring for
    why this hopping pattern gives coordination-1 boundary atoms rather than
    an ordinary zigzag edge's coordination-2).

    Units: energy in units of hopping |t|=1, lengths in nm.

    Reproducibility: this script passes an explicit seed_h to
    kite.Configuration() (see main()'s docstring) so the vacancy placement
    is fixed rather than drawn fresh from the OS's random source every run
    -- re-running it puts the defects at the same sites again. See
    docs/api/kite.md's Configuration entry for what seed_h/seed_v control.

    Last updated: 27/07/2026
"""

import sys
import os

import kite
import numpy as np
from kite import lattice as latt

__all__ = ["afm_zigzag_lattice", "register_sz_operators", "main"]

T = 1.0
DELTA = 0.3

A = 0.24595
A_CC = 0.142
A1 = A * np.array([1.0, 0.0])
A2 = A * np.array([0.5, 0.5 * np.sqrt(3.0)])
POS_A = np.array([0.0, -A_CC / 2])
POS_B = np.array([0.0, A_CC / 2])


def afm_zigzag_lattice(t=T, delta=DELTA):
    """4-sublattice (Aup, Bup, Adn, Bdn) honeycomb lattice with a fixed
    staggered onsite mass +delta on (A, spin up)/(B, spin down),
    -delta on (B, spin up)/(A, spin down) -- i.e. opposite sign for the two
    spins, same alternating-sublattice pattern within each spin."""
    lat = latt.Lattice(a1=A1, a2=A2)
    lat.add_sublattices(
        ('Aup', list(POS_A), +delta),
        ('Bup', list(POS_B), -delta),
        ('Adn', list(POS_A), -delta),
        ('Bdn', list(POS_B), +delta),
    )
    lat.add_hoppings(
        ([0, 0], 'Aup', 'Bup', -t), ([1, -1], 'Aup', 'Bup', -t), ([0, -1], 'Aup', 'Bup', -t),
        ([0, 0], 'Adn', 'Bdn', -t), ([1, -1], 'Adn', 'Bdn', -t), ([0, -1], 'Adn', 'Bdn', -t),
    )
    return lat


def register_sz_operators(calculation):
    """Sz resolved per sublattice (l0=Sz_A, l1=Sz_B) -- a single Sz spanning
    all 4 orbitals would collapse the two-atom basis into one number per
    unit cell and hide the sublattice-alternating (checkerboard) Neel order
    and the edge-resolved sublattice polarization this figure needs."""
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)
    calculation.add_orbital_coupling('Aup', 'Aup', 0.5, 'l0')
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bup', 0.5, 'l1')
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l1')
    return ['l0', 'l1']


def main(t=T, delta=DELTA, lx=8, ly=24, vacancy_concentration=0.0,
         energy=0.05, sigma=0.1, vectors=4000,
         seed_h=987654321,
         output_file=None):
    """vacancy_concentration=0.0 -> clean ribbon. A nonzero value (e.g. 0.05
    for 5%) adds real vacancies: both spin channels removed at the same
    site, since add_vacancy('Aup') and add_vacancy('Adn') on one shared
    StructuralDisorder instance use the same random site selection (one
    random site index per group, both orbitals at that site removed
    together -- see Src/Hamiltonian/HamiltonianVacancies.cpp).

    seed_h fixes the vacancy placement, so re-running this script puts the
    defects at the same sites instead of a fresh random draw each time --
    this is what matters for a clean-vs-vacancy comparison to be about the
    same disorder realization. seed_v (the ldos_map stochastic vectors) is
    left at KITE's default (0, "don't care"): with vectors=4000 the
    stochastic estimate is already self-averaged, so which particular random
    vectors it used doesn't change the physical picture the way which sites
    got vacancies does. Pass seed_h=0 to go back to "don't care" for the
    vacancy positions too."""
    lattice = afm_zigzag_lattice(t, delta)

    struc_disorder = None
    if vacancy_concentration > 0.0:
        struc_disorder = kite.StructuralDisorder(lattice, concentration=vacancy_concentration)
        struc_disorder.add_vacancy('Aup')
        struc_disorder.add_vacancy('Adn')

    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'open'],
        is_complex=True, precision=1, spectrum_range=[-4, 4],
        seed_h=seed_h,
    )
    calculation = kite.Calculation(configuration)
    operators = register_sz_operators(calculation)
    calculation.ldos_map(energy_=energy, sigma_=sigma, vectors_=vectors, operators=operators)

    if output_file is None:
        tag = "clean" if vacancy_concentration == 0.0 else f"vac{int(round(vacancy_concentration * 100))}"
        output_file = f"afm_zigzag_ribbon_{tag}-output.h5"

    kite.config_system(lattice, configuration, calculation,
                        disorder_structural=struc_disorder, filename=output_file)
    return output_file


if __name__ == "__main__":
    conc = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
    fname = main(vacancy_concentration=conc)
    print(f"Wrote {fname}")
    print("Run:   ../build/KITEx", fname)
    print("Then:  python afm_zigzag_ribbon_process.py")
