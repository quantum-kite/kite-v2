""" Post-processing for rashba_edelstein_graphene.py: Kubo-Bastin Edelstein integral.

    This script reads the raw double-Chebyshev moment matrix
    /Calculation/CustomTwo/Gamma written by KITEx, reconstructs the
    Kubo-Bastin integrand using the general (NumVelocities mod 4) formula
    derived in edelstein()'s own docstring, and integrates it against the
    Fermi function to give chi_yx(mu).

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

    General Kubo-Bastin reconstruction for custom_two(), valid for ANY
    NumVelocities parity -- not the spin-Hall-specific ".imag branch" formula
    this file previously copied verbatim from kane_mele_spin_hall_process.py.

    Derivation (verified against both a p=2 case, spin Hall, and this file's
    own p=1 case, pure Rashba):

    Write each raw stored operator as X = i^n_X * X_H with X_H genuinely
    Hermitian (n_X counts "missing factors of i" from KITE's convention that a
    bare commutator/velocity token is stored anti-Hermitian, not with its
    compensating 1/i -- see custom_one's docstring). The physical
    (Hermitian-operator) Kubo-Bastin integrand is
        X(E) = -2*Im[ Tr[A_H delta(E-H) B_H dG+/dE] ],
    and custom_two's stochastic trace gives Gamma_mn = Tr[T_m(H) A T_n(H) B]
    = i^p * Tr[T_m(H) A_H T_n(H) B_H], p = n_A + n_B = NumVelocities. Since
    i^p has period 4 (not 2), the physically correct component of
        Z(E) = sum_mn delta_m(E) * Gamma_mn * dgreen_n(E)
    depends on p mod 4, not just its parity:
        p%4 == 0: X(E) = -2*Im[Z]
        p%4 == 1: X(E) = +2*Re[Z]      (this vertex: bare s_x has 0
                                         velocity tokens, v_y has 1 -> p=1)
        p%4 == 2: X(E) = +2*Im[Z]      (spin Hall: {v_x,s_z}/2 has 1, v_y
                                         has 1 -> p=2 -- this is the case
                                         the old ".imag" formula silently
                                         assumed for every vertex pair)
        p%4 == 3: X(E) = -2*Re[Z]

    Verified: for p=2 (spin Hall), this reduces to exactly (not just
    proportionally) the already-validated kane_mele_spin_hall_process.py
    plateau sigma_xy~-2.02. For p=1 (this vertex), it fixes two independent,
    previously-broken structural checks: chi_yx(mu) becomes odd in mu with a
    genuine zero at mu=0 for pure Rashba (lambda_I=0), and drops to
    near-zero inside the lattice's true bulk gap (lambda_I!=0) instead of
    showing a spurious flat Fermi-sea-like plateau there.

    Units: same e^2/h calibration as conductivity_dc (see custom_two's own
    docstring in docs/api/kite.md), with the EnergyScale^(NumVelocities-2)
    correction applied for NumVelocities != 2 (derived from KITE's rescaled-
    Hamiltonian convention, same logic as custom_one's EnergyScale factor).
    For this vertex (NumVelocities=1) the result carries one fewer power of
    e than e^2/h, i.e. e/h.
    """
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()
        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        energy_scale = np.array(f["EnergyScale"]).item()
        num_velocities = int(np.array(f["/Calculation/CustomTwo/NumVelocities"]))
        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    Moments_D, Moments_G = moments_matrix.shape

    beta = energy_scale / k_BT
    scat_dim = scat_phys / energy_scale
    deltascat_dim = deltascat_phys / energy_scale

    E_grid = np.linspace(-0.995, 0.995, n_egrid)

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    Z = np.einsum("ni,nm,im->i", delta, moments_matrix, dgreenR)

    p = num_velocities % 4
    X = {0: -2.0 * Z.imag, 1: 2.0 * Z.real, 2: 2.0 * Z.imag, 3: -2.0 * Z.real}[p]

    chi = np.empty(len(mu_values))
    for i, mu in enumerate(mu_values):
        integrand = X * fermi_function(E_grid, mu / energy_scale, beta)
        chi[i] = simpson(integrand, E_grid)

    units = 1.0 / (2.0 * np.pi)
    density_scale = num_orbitals / (unit_cell_area * units)
    energy_scale_correction = energy_scale ** (num_velocities - 2)
    chi_yx = 2.0 * chi * density_scale * energy_scale_correction
    return mu_values, chi_yx


def plot(pos_h5, neg_h5, control_h5, out_path="plots/rashba_edelstein_graphene_preview.png"):
    """Compare +lambda_R, -lambda_R, and a lambda_R=0 (Kane-Mele-only) control.

    Flipping the sign of lambda_R must flip the sign of chi_yx (S_REE ~
    zhat x E, and reversing the Rashba term reverses the helicity of the
    spin-momentum locking) -- the same sign-inversion check shown in Fig. 2(a)
    of Medina Duenas et al., Commun. Phys. (2026) [arXiv:2510.21240], just
    with two lambda_R values instead of their full color-graded sweep.

    The control isn't expected to be exactly zero at finite num_random_ (see
    rashba_edelstein_graphene.py's docstring), but should be much smaller than
    -- and structurally distinct from -- the actual signal.
    """
    mus = np.linspace(-2.0, 2.0, 201)
    mu, chi_pos = edelstein(pos_h5, mus)
    _, chi_neg = edelstein(neg_h5, mus)
    _, chi_control = edelstein(control_h5, mus)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(mu, chi_pos, color=kite_style.KITE_PRIMARY, lw=1.5,
            label=r"$\lambda_R=+0.1t$")
    ax.plot(mu, chi_neg, color="#357EDD", lw=1.5,
            label=r"$\lambda_R=-0.1t$")
    ax.plot(mu, chi_control, color=kite_style.KITE_ACCENT, lw=1.5, ls="--",
            label=r"$\lambda_R=0$ (Kane-Mele only, control)")
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
    return mu, chi_pos, chi_neg, chi_control


def plot_pure_rashba(pos_h5, neg_h5,
                      out_path="plots/rashba_edelstein_graphene_pureR_preview.png"):
    """Pure Rashba (lambda_I=0): chi_yx(mu) is odd, with a genuine zero at mu=0.

    No Kane-Mele term means no bulk gap and no band-edge resonance to
    complicate the picture (unlike the full lambda_I!=0 case in plot()) --
    this is the cleanest structural check of the two, and the mu range is
    narrow (-0.5t to 0.5t) since that's where the antisymmetric structure
    actually lives; nothing interesting happens further out for this vertex
    pair without Kane-Mele.
    """
    mus = np.linspace(-0.5, 0.5, 101)
    mu, chi_pos = edelstein(pos_h5, mus)
    _, chi_neg = edelstein(neg_h5, mus)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(mu, chi_pos, color=kite_style.KITE_PRIMARY, lw=1.5,
            label=r"$\lambda_R=+0.1t$")
    ax.plot(mu, chi_neg, color="#357EDD", lw=1.5,
            label=r"$\lambda_R=-0.1t$")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.axvline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel(r"$\mu\,/\,t$", fontsize=13)
    ax.set_ylabel(r"$\chi_{yx}(\mu)$  [$e/h$]", fontsize=13)
    ax.set_title(r"Pure Rashba ($\lambda_I=0$): $\chi_{yx}(\mu)$ is odd, zero at $\mu=0$",
                 fontsize=12)
    ax.legend(loc="best", fontsize=9, framealpha=0.85)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    return mu, chi_pos, chi_neg


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "rashba_edelstein_graphene-output.h5"
    mus = np.linspace(-2.0, 2.0, 201)
    mu, chi_yx = edelstein(fname, mus)
    for m, c in zip(mu, chi_yx):
        print(f"{m: .4f}  {c: .6e}")
