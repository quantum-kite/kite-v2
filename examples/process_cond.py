# Sums a matrix of chebyshev moments, as done in KITE-tools
import h5py
import numpy as np
from scipy.integrate import simpson


def gauss_first(n, mu, sigma):
    n_f = float(n)
    tmp = 1.0 - mu**2
    numerator = n_f * mu * sigma**2 * (n_f**2 * sigma**2 / tmp - 3.0)
    denominator = tmp ** (-1.5)
    return numerator * denominator


def gauss_second(n, mu, sigma):
    n_f = float(n)
    term_1 = 7.0 * mu**2 - 4.0
    tmp = 1.0 - mu**2
    term_2 = 3.0 - 6.0 * n_f**2 * sigma**2 / tmp + (n_f * sigma) ** 4 / (tmp**2)
    denominator = 1.0 / (24.0 * tmp)
    return sigma**2 * term_1 * term_2 * denominator


def build_gaussian(energy_phys, sigma_phys, energy_scale, Moments_G):
    sigma = sigma_phys / energy_scale
    E = energy_phys / energy_scale
    coefs = np.zeros((Moments_G, len(E)), dtype=float)
    coefs[0, :] = 1.0 - gauss_second(0.0, E, sigma)
    for n in range(1, Moments_G):
        n_f = float(n)
        gaussian = np.exp(-0.5 * n_f * n_f * sigma**2 / (1.0 - E**2))
        cossine = np.cos(n_f * np.arccos(E))
        sine = np.sin(n_f * np.arccos(E))

        g1 = gauss_first(n_f, E, sigma)
        g2 = gauss_second(n_f, E, sigma)
        coefs[n, :] = 2.0 * gaussian * (cossine * (1.0 - g2) - 0.5 * sine * g1)
    prefactor = 1.0 / (np.pi * np.sqrt(1.0 - E**2))
    coefs *= prefactor
    return coefs / energy_scale


def build_dgreen(z_phys, energy_scale, Moments_D):
    z = z_phys / energy_scale
    coefs = np.zeros((Moments_D, len(z)), dtype=complex)
    sq = 1.0 - z**2
    sqr = np.sqrt(sq, dtype=complex)
    diff = z - 1j * sqr
    current_power = np.ones_like(z, dtype=complex)
    coefs[0, :] = -1j * z * current_power / (sq * sqr)
    current_power *= diff
    for n in range(1, Moments_D):
        coefs[n, :] = 2.0 * (float(n) - 1j * z / sqr) * current_power / sq
        current_power *= diff
    coefs /= energy_scale * energy_scale
    return coefs


def fermi_function(E, mu, beta):
    return 1.0 / (1.0 + np.exp(np.clip(beta * (E - mu), -100, 100)))


def calculate_conductivity(
    file_path, mu_values, k_BT, E_grid, eta, sigma
):
    mem = 16
    with h5py.File(file_path, "r") as f:
        num_orbitals = np.array(f["NOrbitals"]).item()
        latt_vecs = f["LattVectors"][:]
        unit_cell_area = np.abs(
            latt_vecs[0, 0] * latt_vecs[1, 1] - latt_vecs[0, 1] * latt_vecs[1, 0]
        )
        spin_degeneracy = 1.0
        energy_scale = np.array(f["EnergyScale"]).item()
        moments_matrix = f["/Calculation/CustomTwo/Gamma"][:].T

    pol_g = int(mem * (1 + (18.0 / (eta / energy_scale)) // mem))
    pol_d = int(mem * (1 + (6.0 / (sigma / energy_scale)) // mem))

    delta = build_gaussian(E_grid, sigma, energy_scale, pol_d)
    dgreenR = build_dgreen(E_grid + 1j * eta, energy_scale, pol_g)
    GammaE = np.einsum("mi,mn,ni->i", delta, moments_matrix[:pol_d, :pol_g], dgreenR)
    return GammaE.imag

    # cond_dc = np.zeros(len(mu_values), dtype=complex)
    # beta = energy_scale / k_BT
    # for i, mu in enumerate(mu_values):
    #     mu_scaled = mu / energy_scale
    #     integrand = GammaE * fermi_function(E_grid, mu_scaled, beta)
    #     cond_dc[i] = simpson(integrand, E_grid)

    # units = 1.0 / (2.0 * np.pi)
    # density_scale = (num_orbitals * spin_degeneracy) / (unit_cell_area * units)
    # return 2.0 * cond_dc.imag * density_scale
