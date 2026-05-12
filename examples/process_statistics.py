# import process_cond as cond
# import h5py
# import numpy as np
# import glob as gl
# import sys

# size = 1024
# NEnergies = 512 + 1
# E_min, E_max = -2.0, 2.0
# E_grid = np.linspace(E_min, E_max, NEnergies)
# # E_grid = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
# mu_values = np.array([0.0])
# beta = 0.05

# # eta = 0.8
# # sigma = 0.1
# eta = 0.2
# sigma = 0.05

# prv = np.zeros(len(E_grid))
# avr = np.zeros(len(E_grid))
# var = np.zeros(len(E_grid))
# results = np.zeros((len(E_grid), 2))

# seeds = np.arange(1, 256 + 1)
# count = 0
# for i, seed in enumerate(seeds):
#     all_present = True
#     for j in range(2):
#         path = f"Data/FS{j:02d}_Seed{seed:03d}.h5"
#         with h5py.File(path, "r") as f:
#             item = f.get("/Calculation/CustomTwo/Gamma")
#             if item is not None:
#                 continue
#             else:
#                 all_present = False
#     if all_present:
#         new = np.zeros(len(E_grid))
#         for j in range(2):
#             path = f"Data/FS{j:02d}_Seed{seed:03d}.h5"
#             sign = (-1.0) ** (j % 2)
#             new += sign * cond.calculate_conductivity(path, mu_values, beta, E_grid, eta, sigma)
#         prv = avr.copy()
#         avr += (new - avr) / (count + 1)
#         var += ((new - prv) * (new - avr) - var) / (count + 1)
#         count += 1
# results[:, 0] = avr
# results[:, 1] = np.sqrt(var / count)

# output_path = f"FS_Eta{eta:.3f}_Sigma{sigma:.3f}_Av{count:03d}.dat"
# data_to_save = np.column_stack((E_grid, results))
# np.savetxt(output_path, data_to_save, fmt="%.5e")

import process_cond as cond
import h5py
import numpy as np
import glob as gl
import sys

size = 1024
NEnergies = 512 + 1
E_min, E_max = -2.0, 2.0
# E_grid = np.linspace(E_min, E_max, NEnergies)
E_grid = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
mu_values = np.array([0.0])
beta = 0.05

# eta = 0.8
# sigma = 0.1
eta = 0.2
sigma = 0.05

# seeds = np.arange(1, 32 + 1)
seeds = np.array([1])
avr = np.zeros(len(E_grid))
results = np.zeros((len(E_grid), len(seeds)))

count = 0
for i, seed in enumerate(seeds):
    all_present = True
    for j in range(1):
        # path = f"Data/FS{j:02d}_Seed{seed:03d}.h5"
        path = f"fs.h5"
        with h5py.File(path, "r") as f:
            item = f.get("/Calculation/CustomTwo/Gamma")
            if item is not None:
                continue
            else:
                all_present = False
    if all_present:
        new = np.zeros(len(E_grid))
        for j in range(1):
            path = f"fs.h5"
            sign = (-1.0) ** (j % 2)
            new += sign * cond.calculate_conductivity(path, mu_values, beta, E_grid, eta, sigma)
        results[:, count] = new
        count += 1
print(results)

# output_path = f"FS0_CleanOnes_Eta{eta:.3f}_Sigma{sigma:.3f}.dat"
# data_to_save = np.column_stack((E_grid, results))
# np.savetxt(output_path, data_to_save, fmt="%.5e")
