import cond_sum as cond
import h5py
import numpy as np
import glob as gl
import sys

size = 32
NEnergies = 512 + 1
E_min, E_max = -0.99, 0.99
# E_grid = np.linspace(E_min, E_max, NEnergies)
E_grid = np.array([0.25, 0.0, -0.50]) / 5

mu_values = np.linspace(-4.0, 4.0, 400)
k_BT = 0.05
pol = 32
scat_phys = 0.2
deltascat_phys = 1.0 * scat_phys

results = np.zeros(len(mu_values))
path = f"haldane-output.h5"
new = cond.calculate_conductivity(
    path, mu_values, k_BT, E_grid, scat_phys, deltascat_phys, 1
)
output_path = f"result0.dat"
data_to_save = np.column_stack((mu_values, new))
np.savetxt(output_path, data_to_save, fmt="%.5e")
