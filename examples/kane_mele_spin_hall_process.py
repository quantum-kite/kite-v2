""" Post-processing for kane_mele_spin_hall.py: Kubo-Bastin spin Hall integral.

    There is NO KITE-tools step for custom_two. This script reads the raw
    double-Chebyshev moment matrix /Calculation/CustomTwo/Gamma written by
    KITEx, reconstructs the Kubo-Bastin integrand at a finite reconstruction
    broadening, and integrates it against the Fermi function to give the spin
    Hall conductivity sigma^{s_z}_{xy} as a function of Fermi energy.

    This is the diag=0 (antisymmetric / Hall) branch of examples/cond_sum.py,
    reproduced here so the example is self-contained. Response is returned in
    the same convention as KITE's charge-conductivity post-processors (density
    scaling included); for a spin current the natural quantum is e/(4*pi) but
    the point of this example is the FLAT PLATEAU across the bulk gap and its
    sign, not an exactly quantized number at L=64.

    Usage:
        python kane_mele_spin_hall_process.py kane_mele_spin_hall-output.h5
"""

import sys
import h5py
import numpy as np
from scipy.integrate import simpson


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


def spin_hall(file_path, mu_values, k_BT=0.01, scat_phys=0.04,
              deltascat_phys=0.04, n_egrid=2000):
    """Return (mu_values, sigma^{sz}_xy) using the diag=0 Hall branch."""
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()
        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        energy_scale = np.array(f["EnergyScale"]).item()
        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    Moments_D, Moments_G = moments_matrix.shape
    spin_degeneracy = 1.0

    beta = energy_scale / k_BT
    scat_dim = scat_phys / energy_scale
    deltascat_dim = deltascat_phys / energy_scale

    # Integration grid in normalized energy (dimensionless, -1..1).
    E_grid = np.linspace(-0.995, 0.995, n_egrid)

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    GammaE = np.einsum("ni,nm,im->i", delta, moments_matrix.imag, dgreenR)

    cond = np.zeros(len(mu_values), dtype=complex)
    for i, mu in enumerate(mu_values):
        integrand = GammaE * fermi_function(E_grid, mu / energy_scale, beta)
        cond[i] = simpson(integrand, E_grid)

    units = 1.0 / (2.0 * np.pi)
    density_scale = (num_orbitals * spin_degeneracy) / (unit_cell_area * units)
    return mu_values, 4.0 * cond.real * density_scale


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "kane_mele_spin_hall-output.h5"
    mus = np.linspace(-2.0, 2.0, 201)
    mu, sxy = spin_hall(fname, mus)
    for m, s in zip(mu, sxy):
        print(f"{m: .4f}  {s: .6e}")
