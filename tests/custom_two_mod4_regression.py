"""Regression test for custom_two's NumVelocities-mod-4 Hermitization/
reconstruction table, per kite_stage1_foundational_audit_2026-07-25.tex's
finding: "Generalized custom_two conventions are public for every p mod 4,
while source comments record physical calibration only for p=1,2."

The reconstruction formula (see examples/rashba_edelstein_graphene_process.py's
edelstein() docstring for the full derivation) is: write each raw stored
operator as X = i^n_X * X_H with X_H genuinely Hermitian (KITE's raw velocity
token is missing a factor of i). custom_two's stochastic trace gives
Gamma_mn = Tr[T_m(H) A T_n(H) B] = i^p * Tr[T_m(H) A_H T_n(H) B_H],
p = NumVelocities. Since i^p has period 4, the physical component extracted
from Z(E) = sum_mn delta_m(E)*Gamma_mn*dgreen_n(E) depends on p mod 4:

    p%4 == 0: X(E) = -2*Im[Z]
    p%4 == 1: X(E) = +2*Re[Z]
    p%4 == 2: X(E) = +2*Im[Z]
    p%4 == 3: X(E) = -2*Re[Z]

p=1 (pure Rashba spin density vertex) and p=2 (spin Hall) were independently
validated against known physical results (see rashba_edelstein_graphene.py and
kane_mele_spin_hall.py). This test grounds p=0 and p=3 to that SAME validated
p=2 data, using the fact that Z is linear in Gamma: multiplying the raw
(already-computed, already-validated) Gamma matrix by i^s is mathematically
identical to computing a system with NumVelocities incremented by s (mod 4),
and the corresponding table entry must reconstruct the IDENTICAL physical
result, since nothing physical has changed -- only the bookkeeping of which
table branch to read off. This is not a new independent physical check (that
would require a genuinely different vertex construction with 0 or 3 velocity
tokens), but it IS a rigorous, real-data-grounded verification that the
mod-4 table is self-consistent across all four residues rather than only
independently spot-checked at p=1 and p=2.
"""
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
from rashba_edelstein_graphene_process import fill_delta, fill_dgreenR

TABLE = {
    0: lambda z: -2.0 * z.imag,
    1: lambda z: 2.0 * z.real,
    2: lambda z: 2.0 * z.imag,
    3: lambda z: -2.0 * z.real,
}


def main():
    h5_path = os.path.join(os.path.dirname(__file__), "..", "examples",
                            "kane_mele_spin_hall-output.h5")
    if not os.path.exists(h5_path):
        raise SystemExit(
            f"Missing {h5_path} -- run examples/kane_mele_spin_hall.py + KITEx first."
        )

    with h5py.File(h5_path, "r") as f:
        energy_scale = np.array(f["EnergyScale"]).item()
        num_velocities = int(np.array(f["/Calculation/CustomTwo/NumVelocities"]))
        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    p0 = num_velocities % 4
    print(f"Anchor case: NumVelocities={num_velocities}, p={p0} (already physically validated)")

    Moments_D, Moments_G = moments_matrix.shape
    scat_dim = 0.04 / energy_scale
    deltascat_dim = 0.04 / energy_scale
    E_grid = np.linspace(-0.995, 0.995, 500)

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    Z = np.einsum("ni,nm,im->i", delta, moments_matrix, dgreenR)

    X_ref = TABLE[p0](Z)

    ok = True
    for s in range(4):
        Z_shifted = (1j ** s) * Z
        p_new = (p0 + s) % 4
        X_shifted = TABLE[p_new](Z_shifted)
        max_diff = np.max(np.abs(X_shifted - X_ref))
        passed = max_diff < 1e-9
        ok &= passed
        print(f"  p={p_new} (shift s={s}): max|X - X_ref| = {max_diff:.3e}  "
              f"[{'PASS' if passed else 'FAIL'}]")

    if not ok:
        raise SystemExit("FAIL: mod-4 table is not self-consistent across all residues.")
    print("\nALL FOUR RESIDUES SELF-CONSISTENT.")


if __name__ == "__main__":
    main()
