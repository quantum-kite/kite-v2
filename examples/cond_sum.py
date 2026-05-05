# Sums a matrix of chebyshev moments, as done in KITE-tools
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
        g_m = green_coef(m, z)
        greenR[m, :] = g_m.imag * factor
    return greenR


def fill_dgreenR(E_grid, scat_dim, Moments_D):
    z = E_grid + 1j * scat_dim
    dgreenR = np.zeros((len(E_grid), Moments_D), dtype=complex)
    for m in range(Moments_D):
        factor = 2.0 / (2.0 if m == 0 else 1.0)
        d_m = dgreen_coef(m, z)
        dgreenR[:, m] = d_m * factor
    return dgreenR


def fermi_function(E, mu, beta):
    return 1.0 / (1.0 + np.exp(np.clip(beta * (E - mu), -100, 100)))


def calculate_conductivity(
    file_path, mu_values, k_BT, E_grid, scat_phys, deltascat_phys, diag
):
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()

        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        spin_degeneracy = 1.0

        energy_scale = np.array(f["EnergyScale"]).item()

        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    Moments_D, Moments_G = moments_matrix.shape

    beta = energy_scale / k_BT
    scat_dim = scat_phys / energy_scale
    deltascat_dim = deltascat_phys / energy_scale

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    # Approach 1
    # if diag == 0:
    #     GammaE = np.einsum("im,mn,ni->i", dgreenR, moments_matrix.imag, delta)
    # else:
    #     GammaE = np.einsum("im,mn,ni->i", dgreenR, moments_matrix, delta)
    # Approach 2
    if diag == 0:
        GammaE = np.einsum("ni,nm,im->i", delta, moments_matrix.imag, dgreenR)
    else:
        GammaE = np.einsum("ni,nm,im->i", delta, moments_matrix, dgreenR)
    print(GammaE)
    cond_dc = np.zeros(len(mu_values), dtype=complex)
    for i, mu in enumerate(mu_values):
        mu_scaled = mu / energy_scale
        integrand = GammaE * fermi_function(E_grid, mu_scaled, beta)
        cond_dc[i] = simpson(integrand, E_grid)

    units = 1.0 / (2.0 * np.pi)
    density_scale = (num_orbitals * spin_degeneracy) / (unit_cell_area * units)
    if diag == 0:
        return 4.0 * cond_dc.real * density_scale
    else:
        return 2.0 * cond_dc.imag * density_scale


def calculate_integrand(
    file_path, mu_values, k_BT, E_grid, scat_phys, deltascat_phys, diag
):
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()

        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        spin_degeneracy = 1.0

        energy_scale = np.array(f["EnergyScale"]).item()

        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    Moments_D, Moments_G = moments_matrix.shape

    beta = energy_scale / k_BT
    scat_dim = scat_phys / energy_scale
    deltascat_dim = deltascat_phys / energy_scale

    delta = fill_delta(E_grid, deltascat_dim, Moments_G)
    dgreenR = fill_dgreenR(E_grid, scat_dim, Moments_D)
    if diag == 0:
        GammaE = np.einsum("im,mn,ni->i", dgreenR, moments_matrix.imag, delta)
        return GammaE.real
    else:
        GammaE = np.einsum("im,mn,ni->i", dgreenR, moments_matrix, delta)
        return GammaE.imag
