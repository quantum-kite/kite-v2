""" Post-processing for rashba_edelstein_graphene.py: Kubo-Bastin Edelstein integral.

    There is NO KITE-tools step for custom_two. This script reads the raw
    double-Chebyshev moment matrix /Calculation/CustomTwo/Gamma written by
    KITEx, reconstructs the Kubo-Bastin integrand, and integrates it against
    the Fermi function to give chi_yx(mu).

    Uses the .imag branch of Gamma_mn, same as
    kane_mele_spin_hall_process.py's spin_hall() -- confirmed empirically on
    the simpler examples/rashba_edelstein_square.py model that the .real
    branch is identically zero for this vertex pair (bare s_x, v_y) despite
    its odd combined NumVelocities (opposite parity from spin Hall's).

    Usage:
        python rashba_edelstein_graphene_process.py rashba_edelstein_graphene-output.h5
"""

import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

import kite_style

kite_style.apply()


def green_coef(m, z):
    sqrt_1_z2 = np.sqrt(1.0 - z**2, dtype=complex)
    alpha = z - 1j * sqrt_1_z2
    return -1j * (alpha**m) / sqrt_1_z2


def dgreen_coef(m, z):
    sqrt_1_z2 = np.sqrt(1.0 - z**2, dtype=complex)
    alpha = z - 1j * sqrt_1_z2
    return (alpha**m) / (1.0 - z**2) * (m - 1j * z / sqrt_1_z2)


def fill_delta(E_grid, deltascat_dim, Moments_G):
    z = E_grid + 1j * deltascat_dim
    greenR = np.zeros((Moments_G, len(E_grid)), dtype=complex)
    for m in range(Moments_G):
        factor = -2.0 / ((2.0 if m == 0 else 1.0) * np.pi)
        greenR[m, :] = green_coef(m, z).imag * factor
    return greenR


def fill_dgreenR(E_grid, scat_dim, Moments_D):
    z = E_grid + 1j * scat_dim
    dgreenR = np.zeros((len(E_grid), Moments_D), dtype=complex)
    for m in range(Moments_D):
        factor = 2.0 / (2.0 if m == 0 else 1.0)
        dgreenR[:, m] = dgreen_coef(m, z) * factor
    return dgreenR


def fermi_function(E, mu, beta):
    return 1.0 / (1.0 + np.exp(np.clip(beta * (E - mu), -100, 100)))


def edelstein(file_path, mu_values, k_BT=0.01, scat_phys=0.04,
              deltascat_phys=0.04, n_egrid=2000):
    """Return (mu_values, chi_yx), chi_yx in units of e/h.

    KITE's Kubo-Bastin/Gamma2D reconstruction machinery (the 4.0, the 1/(2*pi)
    inside density_scale, and the internal Chebyshev/Green's-function
    normalizations) is calibrated to output e^2/h for a genuine two-velocity
    (NumVelocities=2) correlator -- verified empirically against a Chern
    insulator's quantized Hall plateau (examples/dos_dccond_haldane.py,
    sigma_xy ~= 1.01-1.02 for C=1, ruling out e^2/hbar, e^2/(pi*h), etc.). This
    vertex pair has NumVelocities=1 (bare s_x has zero, v_y has one -- see
    rashba_edelstein_graphene.py's docstring), one power of EnergyScale short
    of conductivity_dc's/spin Hall's NumVelocities=2 case, which is also
    exactly the one velocity leg (the only place a charge e enters the Kubo-
    Bastin current operator) replaced by the dimensionless s_x operator here.
    So the result carries one fewer power of e than e^2/h, i.e. e/h.
    """
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()
        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        energy_scale = np.array(f["EnergyScale"]).item()
        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    Moments_D, Moments_G = moments_matrix.shape

    beta = energy_scale / k_BT
    scat_dim = scat_phys / energy_scale
    deltascat_dim = deltascat_phys / energy_scale

    E_grid = np.linspace(-0.995, 0.995, n_egrid)

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    GammaE = np.einsum("ni,nm,im->i", delta, moments_matrix.imag, dgreenR)

    chi = np.zeros(len(mu_values), dtype=complex)
    for i, mu in enumerate(mu_values):
        integrand = GammaE * fermi_function(E_grid, mu / energy_scale, beta)
        chi[i] = simpson(integrand, E_grid)

    units = 1.0 / (2.0 * np.pi)
    density_scale = num_orbitals / (unit_cell_area * units)
    chi_yx_raw = 4.0 * chi.real * density_scale  # unchanged Kubo-Bastin/Gamma2D machinery
    chi_yx = chi_yx_raw / energy_scale           # NumVelocities=1 correction (vs. =2 for conductivity_dc)
    return mu_values, chi_yx


def plot(rashba_h5, control_h5, out_path="plots/rashba_edelstein_graphene_preview.png"):
    """Compare the lambda_R!=0 REE response against a lambda_R=0 (Kane-Mele-only) control.

    The control isn't expected to be exactly zero at finite num_random_ (see
    rashba_edelstein_graphene.py's docstring), but should be much smaller than
    -- and structurally distinct from -- the actual signal.
    """
    mus = np.linspace(-2.0, 2.0, 201)
    mu, chi = edelstein(rashba_h5, mus)
    _, chi_control = edelstein(control_h5, mus)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(mu, chi, color=kite_style.KITE_PRIMARY, lw=1.5,
            label=r"$\lambda_R\neq0,\ \lambda_I\neq0$ (Rashba-Edelstein)")
    ax.plot(mu, chi_control, color=kite_style.KITE_ACCENT, lw=1.5, ls="--",
            label=r"$\lambda_R=0,\ \lambda_I\neq0$ (Kane-Mele only, control)")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel(r"$\mu\,/\,t$", fontsize=13)
    ax.set_ylabel(r"$\chi_{yx}(\mu)$  [$e/h$]", fontsize=13)
    ax.set_title("Graphene Rashba-Edelstein effect: spin-density response to $E_y$",
                 fontsize=12)
    ax.legend(loc="best", fontsize=9, framealpha=0.85)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    return mu, chi, chi_control


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "rashba_edelstein_graphene-output.h5"
    mus = np.linspace(-2.0, 2.0, 201)
    mu, chi_yx = edelstein(fname, mus)
    for m, c in zip(mu, chi_yx):
        print(f"{m: .4f}  {c: .6e}")
