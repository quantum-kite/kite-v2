"""Regression matrix for calculation.ldos()/ldos_map()'s operators= feature,
per kite_stage1_foundational_audit_2026-07-25.tex Section 2.4 ("Required
regression matrix"). Five operators on the same exactly-solvable system
(independent two-level cells, H = B*sigma_y, no inter-cell hopping):

  I        : operator result must equal the sum of the two ordinary
             (per-orbital) LDOS values.
  sigma_z  : must equal the DIFFERENCE of the two per-orbital LDOS values
             (sigma_z = diag(+1,-1) = P_up - P_dn).
  sigma_x  : real off-diagonal operator -- invariant under the (b,a)<->(a,b)
             transpose bug (sigma_x^T = sigma_x), so this alone would NOT
             have caught it; included per the audit's table for completeness.
  sigma_y  : imaginary off-diagonal, antisymmetric (sigma_y^T = -sigma_y) --
             THE case the transpose bug flips the sign of. Exact and
             stochastic must agree in sign (see sigma_y_regression.py for the
             standalone, more detailed version of this check).
  P_up     : positive projector diag(1,0) -- must exactly reproduce the plain
             per-orbital LDOS at the 'u' orbital (the identity this feature
             must reduce to when O is a single-orbital projector).

All five run through BOTH the exact (calculation.ldos) and stochastic
(calculation.ldos_map) code paths, in one KITE configuration file.
"""
import glob
import os
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import kite
from kite import lattice as latt

KITEX = os.path.join(os.path.dirname(__file__), "..", "build", "KITEx")
KITE_TOOLS = os.path.join(os.path.dirname(__file__), "..", "build", "KITE-tools")

B = 0.3
OUTPUT_FILE = "operator_ldos_regression-output.h5"

# label -> (2x2 matrix, description)
OPERATORS = {
    'l0': (np.eye(2), "I"),
    'l1': (np.array([[1, 0], [0, -1]]), "sigma_z"),
    'l2': (np.array([[0, 1], [1, 0]]), "sigma_x"),
    'l3': (np.array([[0, -1j], [1j, 0]]), "sigma_y"),
    'l4': (np.array([[1, 0], [0, 0]]), "P_up"),
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

    labels = list(OPERATORS.keys())
    calculation.ldos(energy=[B], num_moments=512, position=[[0, 0], [0, 0]],
                      sublattice=['u', 'd'], num_disorder=1, operators=labels)
    calculation.ldos_map(energy_=B, sigma_=0.05, vectors_=3000, operators=labels)

    kite.config_system(lat, configuration, calculation, filename=OUTPUT_FILE)


def main():
    build()
    subprocess.run([KITEX, OUTPUT_FILE], check=True, capture_output=True)
    subprocess.run([KITE_TOOLS, OUTPUT_FILE, "--LDOS", "-K", "jackson"],
                    check=True, capture_output=True)

    # Plain per-orbital LDOS (position 0 = 'u', position 1 = 'd').
    plain_dat = np.atleast_2d(np.loadtxt(sorted(glob.glob("ldos0.3*.dat"))[0]))
    ldos_u, ldos_d = plain_dat[0, -1], plain_dat[1, -1]

    exact = {}
    for label in OPERATORS:
        files = sorted(glob.glob(f"ldos_{label}_0.3*.dat"))
        exact[label] = np.atleast_2d(np.loadtxt(files[0]))[0, -1]

    with h5py.File(OUTPUT_FILE, "r") as f:
        grp = f["Calculation"]["ldos_map"]["Map_Operators"]
        stoch = {label: np.array(grp[label])[0, 0] for label in OPERATORS}

    ok = True

    def check(name, got, expect, tol):
        nonlocal ok
        passed = abs(got - expect) < tol
        ok &= passed
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}: got={got:+.4f} expect={expect:+.4f}")

    print(f"Plain per-orbital LDOS: u={ldos_u:+.4f}, d={ldos_d:+.4f}\n")

    print("I (exact):")
    check("I == LDOS_u + LDOS_d", exact['l0'], ldos_u + ldos_d, 0.05 * abs(ldos_u + ldos_d) + 1e-6)

    print("sigma_z (exact):")
    check("sigma_z == LDOS_u - LDOS_d", exact['l1'], ldos_u - ldos_d, 0.05 * abs(ldos_u - ldos_d) + 1e-6)

    print("P_up (exact):")
    check("P_up == LDOS_u", exact['l4'], ldos_u, 0.05 * abs(ldos_u) + 1e-6)

    print("sigma_y sign agreement (exact vs stochastic):")
    sign_match = np.sign(exact['l3']) == np.sign(stoch['l3'])
    ok &= sign_match
    print(f"  [{'PASS' if sign_match else 'FAIL'}] exact={exact['l3']:+.4f} stochastic={stoch['l3']:+.4f}")

    print("\nsigma_x (exact, real off-diagonal, transpose-invariant so this alone doesn't test the bug):")
    print(f"  exact={exact['l2']:+.4f} stochastic={stoch['l2']:+.4f}")

    print()
    if not ok:
        raise SystemExit("FAIL: one or more operator-LDOS regression checks failed.")
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
