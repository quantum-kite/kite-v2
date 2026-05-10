import sys
sys.path.append("..")

from interfaces import kite
from interfaces import custom
from interfaces import lattice as latt
import numpy as np

def haldane(onsite=(0, 0), t=1):
    a = 0.24595
    a_cc = 0.142
    t2 = t/10

    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)

    lat.add_sublattices(
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', -t),
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t),
        ([1, 0], 'A', 'A', -t2 * 1j),
        ([0, -1], 'A', 'A', -t2 * 1j),
        ([-1, 1], 'A', 'A', -t2 * 1j),
        ([1, 0], 'B', 'B', -t2 * -1j),
        ([0, -1], 'B', 'B', -t2 * -1j),
        ([-1, 1], 'B', 'B', -t2 * -1j)
    )
    return lat


def main(onsite=(0, 0), t=1):
    """Prepare the input file for KITEx"""
    lattice = haldane(onsite, t)

    nx = ny = 1
    lx = ly = 32
    mode = "open"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-5, 5]
    )
    calculation = kite.Calculation(configuration)
    A = custom.Vertex(32, [[1.0j, "rx"]])
    B = custom.Vertex(64, [[1.0j, "ry"]])
    calculation.custom_two(
        stream_=[A, B],
        num_random_=1,
        num_disorder_=1,
        num_points_=1000,
        temperature_=0.01,
    )

    # configure the *.h5 file
    output_file = "custom.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file

if __name__ == "__main__":
    main()
