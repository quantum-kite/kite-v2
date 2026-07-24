""" Post-processing for spin_precession_simple.py: plot <S(t)>.

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


def load(file_path):
    with h5py.File(file_path, "r") as f:
        timestep = float(np.array(f["Calculation"]["gaussian_wave_packet"]["timestep"]))
        energy_scale = float(np.array(f["EnergyScale"]))
        series = {label: np.array(f["Calculation"]["gaussian_wave_packet"][label]).flatten().real
                  for label in LABELS}
    num_points = len(next(iter(series.values())))
    times = np.arange(num_points) * (timestep / energy_scale)
    return times, series


def plot(file_path, B=0.1, out_path="plots/spin_precession_simple_preview.png"):
    times, series = load(file_path)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for label in ['l0', 'l1', 'l2']:
        ax.plot(times, series[label], lw=1.8, label=f"${LABELS[label]}(t)$, KITE")
    ax.plot(times, 0.5 * np.cos(2 * B * times), lw=1, ls='--', color='k',
            label=r"$\frac{1}{2}\cos(2Bt)$ (exact)")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle S(t)\rangle$")
    ax.set_title("Spin precession in a uniform transverse field")
    ax.legend(loc="best", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    return times, series


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "spin_precession_simple-output.h5"
    times, series = plot(fname)
