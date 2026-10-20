import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# The same graphene Hamiltonian twice: once with the energy origin shifted,
# once without. The moments must reproduce their references in both cases, and
# the DOS KITE-tools reconstructs must agree once the energy axis is shifted
# back -- which is what the shift and rescaling are supposed to guarantee.
tol = 1e-8
ok = True

for run, ref in (("config.h5", "configREF.h5"),
                 ("config_sym.h5", "config_symREF.h5")):
    result = subprocess.run([f"SEED=3 ${{KITEX:-../../build/KITEx}} {run}"],
                            capture_output=True, shell=True)
    if result.returncode != 0:
        print("ERROR")
        exit(0)
    res = compare.compare([run, "/Calculation/dos/MU", ref, "/Calculation/dos/MU"])
    if res[0] > tol:
        ok = False

print("OK" if ok else "Problem")
