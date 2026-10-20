"""Fixture for tests/singleshot_tools_regression.py.

Plain graphene with a deliberately **asymmetric** spectrum range. The asymmetry
is the point: it makes the energy shift nonzero, so a result reported in the
wrong convention is visible. With a symmetric range the shift is zero and any
round-trip error hides completely.

Small and cheap on purpose -- this checks that the energies and conductivities
KITE-tools reports are the ones that were asked for, not that the conductivity
is converged.

Run from this directory:
    python graphene_singleshot_check.py
    ../../build/KITEx graphene_singleshot_check-output.h5
    ../../build/KITE-tools graphene_singleshot_check-output.h5 --SingleShot -X
"""
import kite
import numpy as np
from kite import lattice as latt

__all__ = ["graphene", "ENERGIES", "SPECTRUM_RANGE", "main"]

#: Asymmetric, so that energy_shift != 0. Graphene with t = 1 spans [-3, 3].
SPECTRUM_RANGE = [-3.1, 3.3]

#: The Fermi energies the single-shot solver is asked for, in eV.
ENERGIES = [-1.0, -0.5, 0.0, 0.5, 1.0]


def graphene(t=1.0):
    a = 0.24595
    a_cc = 0.142
    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, 0.5 * np.sqrt(3)])
    lat = latt.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(
        ('A', [0, -a_cc / 2], 0.0),
        ('B', [0, a_cc / 2], 0.0),
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', -t),
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t),
    )
    return lat


def main(t=1.0, moments=256, output_file="graphene_singleshot_check-output.h5"):
    lattice = graphene(t)
    configuration = kite.Configuration(
        divisions=[2, 2],
        length=[64, 64],
        boundaries=["periodic", "periodic"],
        is_complex=False,
        spectrum_range=SPECTRUM_RANGE,
        precision=1,
    )
    calculation = kite.Calculation(configuration)
    calculation.singleshot_conductivity_dc(
        energy=ENERGIES,
        num_moments=moments,
        num_random=1,
        num_disorder=1,
        direction='xx',
        eta=0.1,
    )
    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    print(main())
