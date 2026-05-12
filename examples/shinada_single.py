import sys

sys.path.append("..")

from interfaces import kite
from interfaces import custom
from interfaces import lattice as latt
import numpy as np

def cuprate(t, r, t_prime):
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(
        ("d", [0.0, 0.0], 0.0),
        ("px", [0.5, 0.0], 0.0),
        ("py", [0.0, 0.5], 0.0),
    )
    hop_plus = (t + 1j * r) / 2
    hop_minus = (-t + 1j * r) / 2
    lat.add_hoppings(
        ([0, 0], "d", "px", hop_plus),
        ([-1, 0], "d", "px", hop_minus),
        ([0, 0], "d", "py", hop_plus),
        ([0, -1], "d", "py", hop_minus),
    )
    lat.add_hoppings(
        ([1, 0], "px", "py", -t_prime / 4),
        ([1, -1], "px", "py", t_prime / 4),
        ([0, 0], "px", "py", t_prime / 4),
        ([0, -1], "px", "py", -t_prime / 4),
    )
    return lat


def main():
    t = 1.0
    r = 1.5 * t
    t_nnn = 0.5 * t
    lattice = cuprate(t, r, t_nnn)
    nx = ny = 1
    lx = ly = 32
    W = 0.0 * t
    eta = 0.8 * t
    sigma = 0.1 * t
    energy_scale = 2.7 * t + 0.5 * W
    mem = 16
    pol_A = int(mem * (1 + (18.0 / (eta / energy_scale)) // mem))
    pol_B = int(mem * (1 + (6.0 / (sigma / energy_scale)) // mem))

    mode = "open"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-energy_scale, energy_scale],
    )
    calculation = kite.Calculation(configuration)
    A = custom.Vertex(pol_A, [[1.0j, "vy.rx"], [1.0j, "rx.vy"]])
    B = custom.Vertex(pol_B, [[1.0j, "vy"]])
    calculation.custom_singleshot_two(
        stream_=[A, B],
        num_random_=1,
        num_disorder_=1,
        sigma_=[sigma],
        energies_=[-2.0, -0.5, 0.0, 0.5, 2.0],
        gamma_=[eta],
    )
    output_file = "single.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    main()
