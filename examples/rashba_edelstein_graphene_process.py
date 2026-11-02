""" Post-processing for rashba_edelstein_graphene.py: Kubo-Bastin Edelstein integral.

    This script reads the raw double-Chebyshev moment matrix
    /Calculation/CustomTwo/Gamma written by KITEx, reconstructs the
    Kubo-Bastin integrand using the general (NumVelocities mod 4) formula
    derived in edelstein()'s own docstring, and integrates it against the
    Fermi function to give chi_yx(mu).

    Usage:
        python rashba_edelstein_graphene_process.py rashba_edelstein_graphene-output.h5
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

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
    """Return (mu_values, chi_yx) by running KITE-tools, in units of e/h.

    The reconstruction lives in KITE-tools (--CustomTwo) and nowhere else. This
    used to re-implement it in Python, which meant every convention -- the
    component selected from the contraction, the energy-scale power, the
    carrier charge each velocity operator carries -- had to be kept correct in
    two places at once. They drifted, so there is one implementation now and
    this is a thin wrapper around it.

    k_BT, scat_phys, deltascat_phys and n_egrid are passed through to the
    binary, so the defaults here are the ones the figure was made with rather
    than the binary's own.
    """
    binary = Path(__file__).resolve().parents[1] / "build" / "KITE-tools"
    if not binary.exists():
        raise FileNotFoundError(
            f"{binary} not found -- build KITE-tools before post-processing")

    source = Path(file_path).resolve()
    mu_values = np.asarray(mu_values, dtype=float)
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        local = work / source.name
        shutil.copy(source, local)
        subprocess.run(
            [str(binary), local.name, "--CustomTwo",
             "-T", repr(k_BT), "-S", repr(scat_phys), "-d", repr(deltascat_phys),
             "-E", str(int(n_egrid)),
             "-F", repr(float(mu_values.min())), repr(float(mu_values.max())),
             str(len(mu_values)),
             "-N", "custom_two.dat"],
            cwd=work, check=True, capture_output=True, text=True)
        data = np.loadtxt(work / "custom_two.dat")

    # The binary lays its own Fermi grid over the requested range; put it back
    # on the caller's points so the two curves of the figure share an axis.
    return mu_values, np.interp(mu_values, data[:, 0], data[:, 1])


def plot(pos_h5, neg_h5, control_h5, out_path="plots/rashba_edelstein_graphene_preview.png"):
    """Compare +lambda_R, -lambda_R, and a lambda_R=0 (Kane-Mele-only) control.

    Flipping the sign of lambda_R must flip the sign of chi_yx (S_REE ~
    zhat x E, and reversing the Rashba term reverses the helicity of the
    spin-momentum locking).

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
