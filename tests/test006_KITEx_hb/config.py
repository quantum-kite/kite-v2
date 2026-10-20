import sys
import numpy as np
from kite import lattice as latt
import kite

# The gapped model has half-bandwidth sqrt((3t)^2 + Delta^2) = 3.4354, which
# 3.1 does not contain: the rescaled spectrum reached 1.108 and the Chebyshev
# recursion grew 1.586 per moment. kite.estimate_bounds gives 3.7789 with its
# default margin.
energy_scale = 3.8

def graphene_initial(onsite=(1.673819, -1.673819)):

    theta = np.pi / 3
    a1 = np.array([2 * np.sin(theta), 0])
    a2 = np.array([np.sin(theta), 1 + np.cos(theta)])

    lat = latt.Lattice( a1=a1, a2=a2)

    lat.add_sublattices(
    ('A', [0, 0], onsite[0]),
    ('B', [0, 1], onsite[1]))

    lat.add_hoppings(
    ([0, 0], 'A', 'B', - 1),
    ([0, -1], 'A', 'B', - 1),
    ([1, -1], 'A', 'B', - 1))

    return lat

nx = ny = 2
Lx=Ly=256

lattice = graphene_initial()

configuration = kite.Configuration(divisions=[nx, ny], length=[Lx, Ly], boundaries=["periodic", "periodic"], is_complex=False, precision=0, spectrum_range=[-energy_scale, energy_scale])
calculation = kite.Calculation(configuration)
calculation.conductivity_optical_nonlinear(num_points=1000, num_moments=32, num_random=1, num_disorder=1, direction='yyy', temperature=0.01)
kite.config_system(lattice, configuration, calculation, filename='config.h5')
