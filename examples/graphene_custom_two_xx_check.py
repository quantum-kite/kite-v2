""" Minimal, fast companion example for
tests/custom_two_same_vertex_sign_regression.py.

Plain graphene (2 sublattices, real nearest-neighbour hopping, no SOC/gap --
the simplest possible lattice), bundled with BOTH KITE's dedicated
conductivity_dc(direction='xx') and the generic custom_two([vx,vx]) machinery
in one KITEx run, so both reconstructions of the SAME longitudinal
conductivity come from the identical Hamiltonian/energy scale. Small and
cheap on purpose -- this exists only to check the two methods agree in SIGN,
not to produce a publication-quality conductivity curve.

Run:
    python graphene_custom_two_xx_check.py
    ../build/KITEx graphene_custom_two_xx_check-output.h5
    python ../tests/custom_two_same_vertex_sign_regression.py
"""
import kite
import numpy as np
from kite import custom
from kite import lattice as latt

__all__ = ["graphene", "main"]


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


def main(t=1.0, moments=128, output_file="graphene_custom_two_xx_check-output.h5"):
    lattice = graphene(t)

    nx = ny = 2
    lx = ly = 32

    configuration = kite.Configuration(
        divisions=[nx, ny], length=[lx, ly],
        boundaries=["periodic", "periodic"], is_complex=True, precision=1,
        spectrum_range=[-4, 4],
    )
    calculation = kite.Calculation(configuration)

    calculation.dos(num_points=500, num_moments=moments, num_random=5, num_disorder=1)
    calculation.conductivity_dc(num_points=500, num_moments=moments, num_random=5,
                                 num_disorder=1, direction="xx", temperature=0.01)

    A = custom.Vertex(moments, [[1.0, "vx"]])
    B = custom.Vertex(moments, [[1.0, "vx"]])
    calculation.custom_two(stream_=[A, B], num_random_=5, num_disorder_=1,
                            num_points_=500, temperature_=0.01)

    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    main()
