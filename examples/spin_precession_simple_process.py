""" Post-processing for spin_precession_simple.py: plot <S(t)> and check it.

    Usage:
        python spin_precession_simple_process.py spin_precession_simple-output.h5
"""

import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()

LABELS = {'l0': 'S_x', 'l1': 'S_y', 'l2': 'S_z'}


def transverse_field(f):
    """Read B out of the Hamiltonian, so changing it in main() cannot go unnoticed.

    B is the only on-site term that connects the two orbitals. KITE stores a
    hopping's target as orbital*3**D + sum_d (n_d + 1)*3**d, so the same-cell
    entry of the 'up' row that points at 'down' sits at a known address.
    """
    dim = len(np.array(f['L']))
    address = 1 * 3 ** dim + sum(3 ** d for d in range(dim))
    targets = np.array(f['Hamiltonian']['d'])[0]
    hoppings = np.array(f['Hamiltonian']['Hoppings'])[0]
    scale = float(np.array(f['EnergyScale']))
    return abs(complex(hoppings[targets == address][0])) * scale


def load(file_path):
    with h5py.File(file_path, "r") as f:
        # One entry per wave-packet calculation, so a length-1 array, not a scalar.
        timestep = float(np.array(
            f["Calculation"]["gaussian_wave_packet"]["timestep_physical"]).ravel()[0])
        series = {label: np.array(f["Calculation"]["gaussian_wave_packet"][label]).flatten().real
                  for label in LABELS}
        field = transverse_field(f)
    num_points = len(next(iter(series.values())))
    return np.arange(num_points) * timestep, series, field


def plot(file_path, out_path="plots/spin_precession_simple_preview.png", tol=1e-10):
    times, series, B = load(file_path)

    exact = 0.5 * np.cos(2 * B * times)
    residual = np.abs(series['l2'] - exact).max()
    if residual > tol:
        raise SystemExit(
            "<Sz(t)> departs from (1/2)cos(2Bt) by {:.2e}, above the {:.0e} "
            "tolerance. The precession is exact for this model at any packet "
            "width, so a residual this large means the run or the time axis is "
            "wrong, not that the physics is approximate.".format(residual, tol))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for label in ['l0', 'l1', 'l2']:
        ax.plot(times, series[label], lw=1.8, label=f"${LABELS[label]}(t)$, KITE")
    ax.plot(times, exact, lw=1, ls='--', color='k',
            label=r"$\frac{1}{2}\cos(2Bt)$ (exact)")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle S(t)\rangle$")
    ax.set_title("Spin precession in a uniform transverse field")
    ax.legend(loc="best", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path} (max residual against the exact curve: {residual:.2e})")
    return times, series


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "spin_precession_simple-output.h5"
    times, series = plot(fname)
