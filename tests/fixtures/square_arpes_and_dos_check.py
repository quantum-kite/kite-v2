"""Fixture for tests/arpes_exclusive_regression.py.

A square lattice carrying **two** calculations, ARPES and DOS, in one KITEx
run. Two is the point: `-X` means "only this quantity", and a file holding a
single calculation cannot tell a working `-X` from one that is ignored.

Small and fast -- nothing here is meant to be converged.

Run from this directory:
    python square_arpes_and_dos_check.py
    ../../build/KITEx square_arpes_and_dos_check-output.h5
"""
import kite
import numpy as np
from kite import lattice as latt

__all__ = ["square", "main"]


def square(t=1.0):
    lat = latt.Lattice(a1=np.array([1.0, 0.0]), a2=np.array([0.0, 1.0]))
    lat.add_sublattices(('A', [0, 0], 0.0))
    lat.add_hoppings(
        ([1, 0], 'A', 'A', -t),
        ([0, 1], 'A', 'A', -t),
    )
    return lat


def main(t=1.0, moments=64, output_file="square_arpes_and_dos_check-output.h5"):
    lattice = square(t)
    configuration = kite.Configuration(
        divisions=[2, 2],
        length=[32, 32],
        boundaries=["periodic", "periodic"],
        is_complex=True,
        precision=1,
        spectrum_range=[-5, 5],
    )
    k_path = np.array([[0.0, 0.0], [np.pi / 2, 0.0], [np.pi, 0.0]])

    calculation = kite.Calculation(configuration)
    calculation.arpes(k_vector=k_path, num_moments=moments, weight=[1],
                      num_disorder=1)
    calculation.dos(num_points=128, num_moments=moments, num_random=1,
                    num_disorder=1)

    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    print(main())
