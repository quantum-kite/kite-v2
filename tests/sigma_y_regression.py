"""Regression test for the operator-weighted LDOS transpose bug (audit
kite_stage1_foundational_audit_2026-07-25.tex, Section "Operator-weighted
LDOS: merge blocker and solution").

System: independent two-level cells (a trivial square lattice, no inter-cell
hopping, so every site is exactly solvable), H = B*sigma_y at each site.
O = sigma_y. Since [H, O] = 0 here, the exact eigenstates of H are also
eigenstates of O with eigenvalues +-1, so Tr[O * Im G(r,r,E)] must have the
SAME SIGN as the plain LDOS at the eigenvalue E = +B (both are dominated by
the O=+1 eigenstate there) -- a real-symmetric operator would not detect an
index-transpose bug (O^T = O for Sz/real projectors), but sigma_y is
antisymmetric and imaginary (sigma_y^T = -sigma_y), so a transposed
contraction flips the sign of rho_O relative to the correct answer.

This reproduces the audit's own end-to-end check: exact (calculation.ldos)
and stochastic (calculation.ldos_map) must AGREE IN SIGN once the transpose
bug (orb_mtx(b,a) -> orb_mtx(a,b) in
Src/Simulation/Custom/SimulationLDOSOperator.cpp) is fixed.
"""
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


def build(output_file):
    B = 0.3
    lat = latt.Lattice(a1=[1.0, 0.0], a2=[0.0, 1.0])
    lat.add_sublattices(('u', [0, 0], 0.0), ('d', [0, 0], 0.0))
    # H = B*sigma_y = [[0,-iB],[iB,0]] on the (u,d) 2-level cell, no hopping
    # between cells at all (independent two-level systems, exactly solvable).
    lat.add_hoppings(([0, 0], 'u', 'd', -1j * B))

    lx = ly = 16
    configuration = kite.Configuration(
        divisions=[1, 1], length=[lx, ly], boundaries=['periodic', 'periodic'],
        is_complex=True, precision=1, spectrum_range=[-1, 1],
    )
    calculation = kite.Calculation(configuration)
    calculation.add_orbital_index('u', 0)
    calculation.add_orbital_index('d', 1)
    # sigma_y = [[0,-i],[i,0]] in the (u,d) basis. add_orbital_coupling(start_,
    # last_, c_, label_) sets matrix[row=last_, col=start_] = c_ -- so
    # matrix[0,1]=-i (u row, d col) needs start_='d', last_='u'; matrix[1,0]=+i
    # needs start_='u', last_='d'.
    calculation.add_orbital_coupling('d', 'u', -1j, 'l0')
    calculation.add_orbital_coupling('u', 'd', 1j, 'l0')

    calculation.ldos(energy=[B], num_moments=512, position=[[0, 0]],
                      sublattice=['u'], num_disorder=1, operators=['l0'])
    calculation.ldos_map(energy_=B, sigma_=0.05, vectors_=2000, operators=['l0'])

    kite.config_system(lat, configuration, calculation, filename=output_file)


def main():
    output_file = "sigma_y_regression-output.h5"
    build(output_file)

    subprocess.run([KITEX, output_file], check=True, capture_output=True)
    subprocess.run(
        [KITE_TOOLS, output_file, "--LDOS", "-K", "jackson"],
        check=True, capture_output=True, cwd=os.path.dirname(output_file) or ".",
    )

    # Exact (deterministic) result: read the reconstructed .dat file.
    import glob
    dat_files = sorted(glob.glob("ldos_l0_*.dat"))
    if not dat_files:
        raise SystemExit("No exact-method .dat output found -- KITE-tools reconstruction failed.")
    exact_val = np.loadtxt(dat_files[0])
    exact_val = np.atleast_2d(exact_val)[0, -1]

    # Stochastic result.
    with h5py.File(output_file, "r") as f:
        stoch = np.array(f["Calculation"]["ldos_map"]["Map_Operators"]["l0"])[0, 0]

    print(f"exact rho_O(E=B)      = {exact_val:+.6f}")
    print(f"stochastic rho_O(E=B) = {stoch:+.6f}")

    if np.sign(exact_val) != np.sign(stoch):
        raise SystemExit(
            "FAIL: exact and stochastic operator-weighted LDOS disagree in "
            "sign for sigma_y -- the transpose bug is back."
        )
    print("PASS: exact and stochastic agree in sign for sigma_y.")


if __name__ == "__main__":
    main()
