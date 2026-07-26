"""Analytic stochastic regression for calculation.ldos_map()'s operators=
feature, closing the gap the stage-1 re-audit found in
operator_ldos_regression.py: that script only checks EXACT results
quantitatively and the STOCHASTIC map only for sigma_y's sign, leaving I,
sigma_z, sigma_x, and P_up's stochastic output unasserted.

System: independent two-level cells (no inter-cell hopping), H = B*sigma_y,
same as sigma_y_regression.py/operator_ldos_regression.py. At target energy
E=B, KITE's Gaussian-map estimator is a sum of two Gaussians centered on the
system's two eigenvalues (+B, from the sigma_y=+1 eigenstate; -B, from
sigma_y=-1), each weighted by that eigenstate's O-expectation value:

    g_sigma(x) = exp(-x^2 / (2*sigma^2)) / (sigma * sqrt(2*pi))
    rho_I       = g_sigma(0) + g_sigma(2B)   (both eigenstates contribute weight 1)
    rho_sigma_y = g_sigma(0) - g_sigma(2B)   (+1 at E=+B, -1 at E=-B)
    rho_sigma_x = rho_sigma_z = 0            (sigma_x, sigma_z have zero diagonal
                                               in the sigma_y eigenbasis)
    rho_P_up    = 0.5 * rho_I                (P_up = (I + sigma_y)/2 in this basis,
                                               so its map is the average of the I
                                               and sigma_y expressions above)

These give a closed-form expectation for all five operators' STOCHASTIC output,
not just a sign check, without needing to match the deterministic method's
different (Jackson-kernel) broadening. Run from an isolated temporary
directory so output filenames can't collide with a concurrent or stale run
from operator_ldos_regression.py (both use the same lattice and label names).
"""
import glob
import math
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import kite
from kite import lattice as latt

KITEX = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "KITEx"))

B = 0.3
SIGMA = 0.05
VECTORS = 20000
OUTPUT_FILE = "operator_ldos_analytic_stochastic-output.h5"

# label -> (2x2 matrix, description)
OPERATORS = {
    'l0': (np.eye(2), "I"),
    'l1': (np.array([[1, 0], [0, -1]]), "sigma_z"),
    'l2': (np.array([[0, 1], [1, 0]]), "sigma_x"),
    'l3': (np.array([[0, -1j], [1j, 0]]), "sigma_y"),
    'l4': (np.array([[1, 0], [0, 0]]), "P_up"),
}


def gaussian(x, sigma):
    return math.exp(-x ** 2 / (2.0 * sigma ** 2)) / (sigma * math.sqrt(2.0 * math.pi))


def analytic_expectations():
    g0 = gaussian(0.0, SIGMA)
    g2b = gaussian(2.0 * B, SIGMA)
    rho_I = g0 + g2b
    return {
        'l0': rho_I,
        'l1': 0.0,
        'l2': 0.0,
        'l3': g0 - g2b,
        'l4': 0.5 * rho_I,
    }


def build():
    lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
    lat.add_sublattices(('u', [0, 0], 0.0), ('d', [0, 0], 0.0))
    lat.add_hoppings(([0, 0], 'u', 'd', -1j * B))

    lx = ly = 16
    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-1, 1],
    )
    calculation = kite.Calculation(configuration)
    calculation.add_orbital_index('u', 0)
    calculation.add_orbital_index('d', 1)

    for label, (mat, _) in OPERATORS.items():
        for i, si in enumerate(['u', 'd']):
            for j, sj in enumerate(['u', 'd']):
                if mat[i, j] != 0:
                    calculation.add_orbital_coupling(sj, si, complex(mat[i, j]), label)

    calculation.ldos_map(energy_=B, sigma_=SIGMA, vectors_=VECTORS, operators=list(OPERATORS.keys()))
    kite.config_system(lat, configuration, calculation, filename=OUTPUT_FILE)


def main():
    with tempfile.TemporaryDirectory(prefix="kite_operator_ldos_analytic_") as tmpdir:
        cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            build()
            subprocess.run([KITEX, OUTPUT_FILE], check=True, capture_output=True)

            expected = analytic_expectations()
            with h5py.File(OUTPUT_FILE, "r") as f:
                grp = f["Calculation"]["ldos_map"]["Map_Operators"]
                stoch = {label: np.array(grp[label])[0, 0] for label in OPERATORS}
                stoch_err = {label: np.array(grp[label])[1, 0] for label in OPERATORS}
        finally:
            os.chdir(cwd)

    print(f"Analytic reference: sigma={SIGMA}, B={B}")
    ok = True

    def check(label):
        nonlocal ok
        name = OPERATORS[label][1]
        got, err, want = stoch[label], stoch_err[label], expected[label]
        # 6-sigma tolerance on the stochastic standard error, floored so that
        # an exactly-zero expectation (sigma_x, sigma_z) isn't checked with a
        # zero-width window.
        tol = max(6.0 * err, 1e-3)
        passed = abs(got - want) < tol
        ok &= passed
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}: got={got:+.5f} +/- {err:.5f}, "
              f"analytic={want:+.5f}")

    for label in OPERATORS:
        check(label)

    print()
    if not ok:
        raise SystemExit("FAIL: one or more analytic stochastic operator-LDOS checks failed.")
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
