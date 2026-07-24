""" Post-processing for oam_quadrupole_precession.py: plot <L(t)> and <Q(t)>.

    Usage:
        python oam_quadrupole_precession_process.py oam_quadrupole_precession-output.h5
"""

import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()

LABELS = {
    'l0': 'L_x', 'l1': 'L_y', 'l2': 'L_z',
    'l3': 'Q_xx', 'l4': 'Q_yy', 'l5': 'Q_zz',
    'l6': 'Q_xy', 'l7': 'Q_xz', 'l8': 'Q_yz',
}


def load(file_path):
    with h5py.File(file_path, "r") as f:
        timestep = float(np.array(f["Calculation"]["gaussian_wave_packet"]["timestep"]))
        energy_scale = float(np.array(f["EnergyScale"]))
        series = {}
        for label in LABELS:
            series[label] = np.array(f["Calculation"]["gaussian_wave_packet"][label]).flatten().real

    num_points = len(next(iter(series.values())))
    deltaT = timestep / energy_scale
    times = np.arange(num_points) * deltaT
    return times, series


def plot(file_path, out_path="plots/oam_quadrupole_precession_preview.png"):
    times, series = load(file_path)

    fig, (ax_l, ax_q) = plt.subplots(1, 2, figsize=(12, 4.5))

    for label in ['l0', 'l1', 'l2']:
        ax_l.plot(times, series[label], lw=1.5, label=f"${LABELS[label]}$")
    ax_l.set_xlabel(r"$t$")
    ax_l.set_ylabel(r"$\langle L(t)\rangle$")
    ax_l.set_title("Orbital angular momentum")
    ax_l.legend(loc="best", fontsize=9)

    for label in ['l3', 'l4', 'l5', 'l6', 'l7', 'l8']:
        ax_q.plot(times, series[label], lw=1.5, label=f"${LABELS[label]}$")
    ax_q.set_xlabel(r"$t$")
    ax_q.set_ylabel(r"$\langle Q(t)\rangle$")
    ax_q.set_title("Orbital quadrupoles")
    ax_q.legend(loc="best", fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    return times, series


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "oam_quadrupole_precession-output.h5"
    times, series = plot(fname)
    print("Q_xy range (should be exactly 0, a symmetry-protected check):",
          series['l6'].min(), series['l6'].max())
