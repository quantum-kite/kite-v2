"""Regression test for a real sign bug found and fixed in
tools/Src/Custom/customtwo.cpp (2026-08): the "Gamma orientation" sign
correction applied for even NumVelocities (p%4 in {0,2}) was calibrated
ONLY against antisymmetric, A!=B "Hall-type" vertex pairs (spin Hall:
A={v_x,s_z}/2, B=v_y -- see custom_two_mod4_regression.py) and silently
assumed to hold for every vertex pair sharing that parity. It does not: for
a SAME-operator pair (A=B=v_x, the ordinary longitudinal conductivity
sigma_xx), the old fixed sign_convention=-1 gave a result anti-correlated
(r~-0.95) with KITE's independently-implemented conductivity_dc(direction=
'xx') on the identical lattice. The fix distinguishes the two cases at
runtime via Gamma's own realness (A=B pairs produce a real, not genuinely
complex, Hermitized Gamma -- see customtwo.cpp's comment for the full
cyclic-trace argument), which does not touch the (already independently
validated) A!=B branch at all.

This test is the permanent guard against regressing that fix: it runs the
same longitudinal conductivity through BOTH the trusted conductivity_dc(xx)
path and the generic custom_two([vx,vx]) path on a minimal, fast lattice,
and asserts the two agree in sign (positive correlation), not just that the
mod-4 table is internally self-consistent (that is
custom_two_mod4_regression.py's job, and does not catch this bug -- it only
checks relative consistency across a single already-computed Gamma, never
compares against an independent reference).

Usage:
    python examples/graphene_custom_two_xx_check.py
    ../build/KITEx examples/graphene_custom_two_xx_check-output.h5
    python tests/custom_two_same_vertex_sign_regression.py
"""
import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(__file__)
H5_PATH = os.path.join(HERE, "..", "examples", "graphene_custom_two_xx_check-output.h5")
KITE_TOOLS = os.path.join(HERE, "..", "build", "KITE-tools")


def _run_reconstruction(flag, extra_args, out_name):
    out_path = os.path.join(HERE, out_name)
    cmd = [KITE_TOOLS, H5_PATH, flag, *extra_args, "-N", out_path, "-X"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"{' '.join(cmd)} failed:\n{result.stdout}\n{result.stderr}")
    return np.loadtxt(out_path)


def main():
    if not os.path.exists(H5_PATH):
        raise SystemExit(
            f"Missing {H5_PATH} -- run examples/graphene_custom_two_xx_check.py "
            f"+ KITEx first."
        )
    if not os.path.exists(KITE_TOOLS):
        raise SystemExit(f"Missing {KITE_TOOLS} -- build KITE-tools first.")

    with h5py.File(H5_PATH, "r") as f:
        num_velocities = int(np.array(f["/Calculation/CustomTwo/NumVelocities"]))
        gamma = f["/Calculation/CustomTwo/Gamma"][()]

    print(f"NumVelocities = {num_velocities} (expect 2, from vx.vx)")
    imag_weight = np.abs(gamma.imag).sum()
    real_weight = np.abs(gamma.real).sum()
    real_fraction = real_weight / (real_weight + imag_weight)
    print(f"Gamma real fraction = {real_fraction:.3f} (expect > 0.5 for an A=B pair)")
    if real_fraction <= 0.5:
        raise SystemExit(
            "FAIL: Gamma is not real-dominated -- the A=B detection this fix relies "
            "on would not trigger; the test lattice/vertex may have changed."
        )

    common = ["-F", "-3", "3", "60", "-E", "1000", "-S", "0.05", "-d", "0.05", "-T", "0.05"]
    conddc = _run_reconstruction("--CondDC", common, "conddc_xx_regression_check.dat")
    customtwo = _run_reconstruction("--CustomTwo", common, "customtwo_xx_regression_check.dat")

    corr = np.corrcoef(conddc[:, 1], customtwo[:, 1])[0, 1]
    print(f"corr(conductivity_dc(xx), custom_two(vx,vx)) = {corr:+.4f} (expect > +0.7)")

    if corr < 0.7:
        raise SystemExit(
            f"FAIL: custom_two(vx,vx) disagrees in sign/shape with the trusted "
            f"conductivity_dc(xx) reference (corr={corr:+.4f}) -- the same-vertex "
            f"sign fix in customtwo.cpp appears to have regressed."
        )
    print("\nPASS: custom_two(vx,vx) matches conductivity_dc(xx) in sign.")


if __name__ == "__main__":
    main()
