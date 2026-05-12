import process_cond as cond
import h5py
import numpy as np
import glob as gl
import sys

size = 1024
NEnergies = 512 + 1
E_grid = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])

W = 1.0
eta = 0.2
sigma = 0.05

base_dir = "../slurm_tools/Data/"
seeds = np.arange(1, 256 + 1)

prv = np.zeros(len(E_grid))
avr = np.zeros(len(E_grid))
var = np.zeros(len(E_grid))
results = np.zeros((len(E_grid), 2))

count = 0
for i, seed in enumerate(seeds):
    all_present = True
    for j in range(2):
        path = base_dir + f"Shinada{j:01d}_Seed{seed:02d}.h5"
        with h5py.File(path, "r") as f:
            item = f.get("/Calculation/CustomSingleTwo/Result")
            if item is not None:
                continue
            else:
                all_present = False
    if all_present:
        # new = np.zeros(len(mu_values))
        new = np.zeros(len(E_grid))
        for j in range(2):
            path = base_dir + f"Shinada{j:01d}_Seed{seed:02d}.h5"
            sign = (-1.0) ** (j % 2)
            with h5py.File(path, "r") as f:
                data = f["/Calculation/CustomSingleTwo/Result"][:].flatten()
            new += sign * data.imag
        prv = avr.copy()
        avr += (new - avr) / (count + 1)
        var += ((new - prv) * (new - avr) - var) / (count + 1)
        count += 1
print(count)
results[:, 0] = avr
results[:, 1] = np.sqrt(var / count)

output_path = f"Singleshot_Eta{eta:.3f}_Sigma{sigma:.3f}.dat"

data_to_save = np.column_stack((E_grid, results))
np.savetxt(output_path, data_to_save, fmt="%.5e")
