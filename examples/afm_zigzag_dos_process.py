""" Post-processing for afm_zigzag_dos.py: spin-resolved total DOS of a
    zigzag-terminated Neel/AFM graphene ribbon, clean vs. a low (~5%)
    concentration of real vacancies.

    ##########################################################################
    #                         Copyright 2026, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Reads the plain-text KITE-tools output (Tr[P_spin * delta(E-H)] via
    --CustomOne, two columns: energy, spectral density) for four
    calculations -- {clean, 5% vacancy} x {spin up, spin down} -- produced
    by running afm_zigzag_dos.py's four generated .h5 configs through KITEx
    then:
        KITE-tools <file>.h5 --DOS -N dos_<tag>.dat \\
                             --CustomOne -E -4 4 1000 -N custom_<tag>.dat

    Design: two panels (clean | vacancy), NOT four/overlaid-in-one -- the
    physically important comparison is spin-up vs. spin-down WITHIN each
    disorder condition (do they coincide, or not), so each panel gets its
    own pair of curves with a shared, identical y-axis scale across both
    panels (mandatory for a fair "did disorder change this" comparison, same
    reasoning as the real-space companion figure). Vertical dotted guides at
    E=+/-Delta mark the flat-edge-band energy where the local real-space
    figure shows the strongest polarization.

    Usage:
        python afm_zigzag_dos_process.py
    (reads custom_clean_up.dat, custom_clean_down.dat, custom_vac5_up.dat,
    custom_vac5_down.dat from the current directory)
"""

import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()

DELTA = 0.3


def load(tag):
    return np.loadtxt(f"custom_{tag}.dat")


def plot_dos_comparison(out_path="plots/afm_zigzag_dos_preview.png"):
    up_clean, dn_clean = load("clean_up"), load("clean_down")
    up_vac, dn_vac = load("vac5_up"), load("vac5_down")

    up_color = kite_style.KITE_PRIMARY
    dn_color = kite_style.KITE_ACCENT
    ymax = 1.08 * max(up_clean[:, 1].max(), dn_clean[:, 1].max(),
                       up_vac[:, 1].max(), dn_vac[:, 1].max())

    fig, (ax_clean, ax_vac) = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)

    for ax, up, dn, title in [
        (ax_clean, up_clean, dn_clean, "Clean ribbon"),
        (ax_vac, up_vac, dn_vac, "5% vacancy conc. (A sublattice)"),
    ]:
        ax.plot(up[:, 0], up[:, 1], color=up_color, lw=1.8, label=r"spin up")
        ax.plot(dn[:, 0], dn[:, 1], color=dn_color, lw=1.8, ls="--", label=r"spin down")
        for x in (+DELTA, -DELTA):
            ax.axvline(x, color="0.6", lw=0.8, ls=":", zorder=0)
        ax.axvline(0, color="0.85", lw=0.7, zorder=0)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(0, ymax)
        ax.set_xlabel(r"$E$ ($t$)")
        ax.set_title(title)
        ax.grid(False)

    ax_clean.text(-DELTA - 0.08, ymax * 0.93, r"$-\Delta$", fontsize=9.5, color="0.4",
                  ha="right", va="center")
    ax_clean.text(DELTA + 0.08, ymax * 0.93, r"$+\Delta$", fontsize=9.5, color="0.4",
                  ha="left", va="center")
    ax_clean.set_ylabel(r"$\mathrm{Tr}\!\left[P_\sigma\,\delta(E-H)\right]$ (states / $t$)")
    ax_clean.legend(loc="upper left", handlelength=1.6, fontsize=10)

    fig.suptitle("Spin-resolved total DOS: global cancellation vs. vacancy-broken symmetry",
                 fontsize=13.5, fontweight="bold", y=1.02)
    fig.text(0.5, 0.925,
              r"Clean ribbon: spin-up $\equiv$ spin-down globally (edge polarization is purely "
              "local — see the real-space figure). Vacancies (A sublattice only) break that "
              "exact cancellation.",
              ha="center", fontsize=9.5, color="0.35")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    plot_dos_comparison()
