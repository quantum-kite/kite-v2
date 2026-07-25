""" Spin-resolved total density of states of a zigzag-terminated graphene
    ribbon with a fixed staggered Neel/AFM mass term, clean vs. a low (~5%)
    concentration of real vacancies.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    See afm_zigzag_bands.py's module docstring for the full physical picture.
    This script computes the GLOBAL counterpart of afm_zigzag_ribbon.py's
    LOCAL real-space probe: Tr[P_up * delta(E-H)] and Tr[P_dn * delta(E-H)],
    the total (whole-ribbon, not site-resolved) spin-up and spin-down DOS,
    via calculation.custom_one() with the spin-up/spin-down PROJECTORS
    (NOT Sz!) registered separately so they can be compared directly rather
    than only through their difference.

    Expected result: in the CLEAN ribbon, spin-up DOS == spin-down DOS at
    every energy (the two edges are individually spin-polarized -- see the
    companion real-space figure -- but related to each other by an exact
    A<->B, up<->down symmetry of the clean Hamiltonian, so their
    contributions cancel in the *global* trace). Vacancies break that exact
    cancellation (they are only introduced on the A sublattice here -- see
    afm_zigzag_ribbon.py's register/main), so the disordered ribbon's
    spin-up and spin-down DOS are expected to visibly split, especially near
    the edge-band energy E=+/-Delta.

    Units: energy in units of hopping |t|=1, lengths in nm.
    Last updated: 25/07/2026
"""

import sys

import kite
from kite import custom

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from afm_zigzag_ribbon import afm_zigzag_lattice, T, DELTA

__all__ = ["register_spin_projectors", "main"]


def register_spin_projector(calculation, spin):
    """Tr[P_spin*delta(E-H)] -- the spin PROJECTOR (weight 1.0 on that spin's
    two sublattices only), not Sz (which would give only their weighted
    difference). KITE-tools' --CustomOne reconstructs a SINGLE registered
    vertex per h5 file (confirmed against this KITE build: a second
    custom_one() call on the same Calculation is silently not reconstructed
    by KITE-tools), so spin-up and spin-down are written to two SEPARATE
    files rather than as two vertices in one -- see main()."""
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)
    if spin == "up":
        calculation.add_orbital_coupling('Aup', 'Aup', 1.0, 'l0')
        calculation.add_orbital_coupling('Bup', 'Bup', 1.0, 'l0')
    elif spin == "down":
        calculation.add_orbital_coupling('Adn', 'Adn', 1.0, 'l0')
        calculation.add_orbital_coupling('Bdn', 'Bdn', 1.0, 'l0')
    else:
        raise ValueError("spin must be 'up' or 'down'")
    return ['l0']


def main(spin, t=T, delta=DELTA, lx=8, ly=24, vacancy_concentration=0.0,
         num_moments=512, num_random=100, num_disorder=None,
         output_file=None):
    lattice = afm_zigzag_lattice(t, delta)

    struc_disorder = None
    if vacancy_concentration > 0.0:
        struc_disorder = kite.StructuralDisorder(lattice, concentration=vacancy_concentration)
        struc_disorder.add_vacancy('Aup')
        struc_disorder.add_vacancy('Adn')
    if num_disorder is None:
        num_disorder = 1 if vacancy_concentration == 0.0 else 20

    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'open'],
        is_complex=True, precision=1, spectrum_range=[-4, 4],
    )
    calculation = kite.Calculation(configuration)
    register_spin_projector(calculation, spin)

    calculation.dos(num_points=1000, num_moments=num_moments,
                     num_random=num_random, num_disorder=num_disorder)

    vertex = custom.Vertex(num_moments, [[1.0, "l0"]])
    calculation.custom_one(stream_=vertex, num_random_=num_random, num_disorder_=num_disorder)

    if output_file is None:
        tag = "clean" if vacancy_concentration == 0.0 else f"vac{int(round(vacancy_concentration * 100))}"
        output_file = f"afm_zigzag_dos_{tag}_{spin}-output.h5"

    kite.config_system(lattice, configuration, calculation,
                        disorder_structural=struc_disorder, filename=output_file)
    return output_file


if __name__ == "__main__":
    spin_arg = sys.argv[1] if len(sys.argv) > 1 else "up"
    conc = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    fname = main(spin_arg, vacancy_concentration=conc)
    print(f"Wrote {fname}")
    print("Run:   ../build/KITEx", fname)
    print("Then:  ../build/KITE-tools", fname,
          "--DOS -N dos.dat --CustomOne -E -4 4 1000 -N custom_one.dat")
