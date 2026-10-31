""" Minimal companion fixture for tests/second_order_full_path_regression.py.

Gapped honeycomb (hexagonal boron nitride), the same model and the same
`xxy` component as examples/hbn_optcond2_vacancies.py, asking for the general
second-order path -- all four Chebyshev objects Gamma0/1/2/3 -- rather than the
internal two-object shortcut. Small and cheap on purpose: this exists to check
that the full path runs and stays reproducible, not to produce a converged
second-order conductivity.

Seeds are fixed so two independent exports give bit-identical output.

Run from this directory:
    python hbn_optcond2_full_path_check.py
    ../../build/KITEx hbn_optcond2_full_path_check-output.h5
"""
import sys

import numpy as np

import kite
from kite import lattice as latt

__all__ = ["hbn", "main"]

#: As in the example: a is the lattice constant, a_cc the B-N distance.
A, A_CC, T, GAP = 0.246, 0.142, -2.8, 4.5


def hbn(a=A, a_cc=A_CC, t=T, gap=GAP):
    lat = latt.Lattice(a1=[a / 2, a / 2 * np.sqrt(3)],
                       a2=[a / 2, -a / 2 * np.sqrt(3)])
    lat.add_sublattices(('A', [0, 0], -gap / 2.0),
                        ('B', [0, a_cc], gap / 2.0))
    lat.add_hoppings(([0, 0], 'A', 'B', t),
                     ([0, 1], 'A', 'B', t),
                     ([-1, 0], 'A', 'B', t))
    return lat


def main(moments=64, lx=64, output_file="hbn_optcond2_full_path_check-output.h5"):
    lattice = hbn()

    configuration = kite.Configuration(
        divisions=[2, 2], length=[lx, lx],
        boundaries=["periodic", "periodic"], is_complex=True, precision=1,
        spectrum_range=[-9.5, 9.5], seed_h=20261031, seed_v=20261031,
    )
    calculation = kite.Calculation(configuration)
    calculation.conductivity_optical_nonlinear(
        "xxy", num_points=64, num_moments=moments, num_random=1, num_disorder=1)

    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    main(output_file=sys.argv[1] if len(sys.argv) > 1
         else "hbn_optcond2_full_path_check-output.h5")
