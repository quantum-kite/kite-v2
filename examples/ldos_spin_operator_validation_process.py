""" Post-processing / cross-validation for ldos_spin_operator_validation.py.

    Compares the exact (Chebyshev-moment) and stochastic (Markov-map)
    operator-weighted local spectral densities against each other, and runs
    the single-orbital-projector regression check against the plain
    (operators=None) LDOS/Map output. See ldos_spin_operator_validation.py's
    module docstring for the full physics background.

    Usage (after running KITEx and KITE-tools --LDOS on the .h5 file this
    script's companion example produces):
        python ldos_spin_operator_validation_process.py ldos_spin_operator_validation-output.h5
"""

import sys
import glob
import h5py
import numpy as np


def load_stochastic(file_path):
    with h5py.File(file_path, "r") as f:
        grp = f["Calculation"]["ldos_map"]
        # HDF5 datasets are written transposed relative to the Eigen arrays
        # that produced them: row 0 = mean over sites, row 1 = stderr.
        plain = np.array(grp["Map"])            # (2, Sizet)
        ops = {}
        if "Operators" in grp:
            labels = [s.decode() if isinstance(s, bytes) else s for s in grp["Operators"][:]]
            for label in labels:
                ops[label] = np.array(grp["Map_Operators"][label])  # (2, NumSites)
    return plain, ops


def load_exact_dat(prefix, label):
    """Load whichever *_<label>_<energy>.dat files KITE-tools wrote."""
    files = sorted(glob.glob(f"{prefix}_{label}_*.dat"))
    data = {}
    for fname in files:
        energy = float(fname[len(f"{prefix}_{label}_"):-len(".dat")])
        data[energy] = np.atleast_2d(np.loadtxt(fname))
    return data


def main(file_path="ldos_spin_operator_validation-output.h5", dat_prefix="ldos", lx=32):
    plain, ops = load_stochastic(file_path)
    n_sites = plain.shape[1] // 4  # 4 sublattices per unit cell; site 0 = Aup

    print("=== Stochastic method: projector regression check (l1 = Aup projector) ===")
    proj_mean, proj_err = ops['l1'][0], ops['l1'][1]
    plain_aup = plain[0, :n_sites]
    diff = np.abs(proj_mean - plain_aup)
    combined_err = np.sqrt(proj_err ** 2 + plain[1, :n_sites] ** 2)
    frac_ok = np.mean(diff < 5 * (combined_err + 1e-12))
    print(f"  mean|diff| = {diff.mean():.4g}, mean reported stderr = {proj_err.mean():.4g}")
    print(f"  fraction of sites within 5-sigma combined error: {frac_ok:.2%}")

    print("\n=== Sz (l0) real-space texture around the magnetic impurity ===")
    exact_l0 = load_exact_dat(dat_prefix, "l0")
    if not exact_l0:
        print(f"  No {dat_prefix}_l0_*.dat files found -- run KITE-tools --LDOS first.")
        return

    stoch_map = ops['l0'][0]   # (NumSites,) mean, indexed by global site = x + y*lx
    stoch_err = ops['l0'][1]

    for energy, data in exact_l0.items():
        print(f"\n  E={energy:+.3f}:")
        print("  dist-from-impurity   exact Sz     stochastic Sz")
        # Rows are in the order requested (increasing distance from the
        # impurity along x); recover (x,y) from the .dat file and look up
        # the SAME site in the full stochastic map for a direct comparison.
        x0, y0 = int(data[0, 0]), int(data[0, 1])
        for row in data:
            x, y = int(row[0]), int(row[1])
            dist = min((x - x0) % lx, (x0 - x) % lx)
            exact_val = row[-1]
            site_idx = x + y * lx
            print(f"    {dist:>3d}              {exact_val:+.5f}     "
                  f"{stoch_map[site_idx]:+.5f} +/- {stoch_err[site_idx]:.5f}")

    print("\nNote: Sz decays away from the impurity -- this is the real-space texture the")
    print("clean altermagnet cannot show (see the module docstring): the impurity's onsite")
    print("shift acts on 'Aup' only, breaking the (rotation)x(spin-swap) symmetry that")
    print("otherwise forces Sz(r)=0 everywhere, so a local spin polarization forms and")
    print("decays away from the defect, Friedel-oscillation style.")


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "ldos_spin_operator_validation-output.h5"
    main(fname)
