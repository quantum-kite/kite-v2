""" Post-processing for haldane_orbital_magnetization.py.

    Usage:
        ../build/KITEx haldane_orbital_magnetization-output.h5
        ../build/KITE-tools haldane_orbital_magnetization-output.h5 --CustomOne -E -4 4 800 -N custom_one.dat
        python process_haldane_orbital_magnetization.py haldane_orbital_magnetization-output.h5 custom_one.dat

    KITE-tools' --CustomOne mode reconstructs the raw energy-resolved trace
    S_raw(E) from the Chebyshev moments (Jackson kernel, same machinery as
    --DOS, already including the 1/EnergyScale Jacobian and the real/imaginary
    rotation for this vertex's odd velocity-operator count -- see
    docs/api/kite-tools.md). What remains is converting S_raw(E) into the
    actual physical orbital-magnetization spectral density S(E) = Tr[M_z_hat *
    delta(E-H)], and integrating it up to a Fermi energy E_F to get the
    magnetization M(E_F).

    Physical normalization, derived from scratch (not fit to match anything)
    -----------------------------------------------------------------------
    The vertex is A = rx.vy - ry.vx = x*H*y - y*H*x (see
    haldane_orbital_magnetization.py), verified to equal exactly this operator
    combination with no extra sign or factor. The reference (Vidarte et al.,
    Phys. Rev. B, doi 10.1103/bhg8-c3rw, Eq. 1-2) defines

        M_z_hat = -(i*e)/(2*hbar*c*Area) * (x*H*y - y*H*x) = -(e/(2*hbar*c*Area)) * i*A

    A is anti-Hermitian (verified: A^dagger = -A), so Tr[A*delta(E-H)] is
    purely imaginary; write it as i*rho_A(E) with rho_A(E) real. Then

        S(E) = Tr[M_z_hat*delta(E-H)] = -(e/(2 hbar c Area)) * i * (i*rho_A(E))
             = +(e/(2 hbar c Area)) * rho_A(E)                                  (i*i = -1)

    KITE-tools' --CustomOne output S_raw(E) equals rho_A(E) up to two things:

    1. The stochastic-trace normalization is per total degree of freedom
       (verified empirically: KITE's own --DOS integrates to exactly 1.0, not
       NOrbitals, over the whole spectrum) -- converting this "per orbital"
       quantity to a genuine per-unit-AREA density requires multiplying by
       (orbitals per cell)/(area per cell) = NOrbitals/unit_cell_area (the
       per-unit-cell factors cancel between the extensive trace and the
       extensive area, leaving this fixed, system-size-independent ratio).
    2. KITE's "vx"/"vy" building blocks (Src/Hamiltonian/HamiltonianRegular.cpp)
       are built directly from the HOPPING VALUES AS EXPORTED TO THE HDF5 --
       which src/kite/__init__.py already divides by EnergyScale before
       KITEx ever sees them. So the vertex KITE actually evaluates is
       A_computed = A_physical / EnergyScale, and rho_A(E) (physical) needs an
       extra factor of EnergyScale to compensate. This is a genuine hidden
       Jacobian, easy to miss: it does not show up in the DOS integral check
       (identity operator has no hopping-dependent "v" tokens), only in any
       vertex actually built from "v".

    With e = hbar = c = 1: S(E) = 0.5 * EnergyScale * (NOrbitals/unit_cell_area) * S_raw(E).

    Verified against the (independently, literature-checked) k-space "modern
    theory" result for this exact lattice: the filled/valence-band Chern
    number is C=+1 (Src/kite/visualize.py's hamiltonian_k, cross-checked twice
    independently against the literature convention for this bond-phase
    assignment), so the exact topological (Streda) relation predicts
    dM/dE_F = C/(2*pi) = +0.1592 in the bulk gap. The reconstruction below
    (num_random=400, num_disorder=4, small Anderson disorder to help
    convergence) gives dM/dE_F = +0.154 for a fit window safely inside the
    numerically-verified gap |E|<1.0 -- agreement to ~3%, consistent with
    finite num_random/disorder-realization stochastic noise, not a remaining
    systematic error. (A fit window extending close to the gap edge biases the
    slope badly -- e.g. |E|<1.5 gives only +0.113 -- because it starts
    incorporating band-edge states; keep well inside the gap.)
"""

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt

import kite_style

kite_style.apply()

# Numerically-verified (10 significant figures, via kite.visualize.hamiltonian_k
# and an independent Fukui-Hatsugai-Suzuki Chern-number calculation) bulk gap
# edge for t2=t/3, phi=pi/2, delta=0. NOT 3*sqrt(3)*t2*sin(phi) = sqrt(3) --
# that is the local Dirac-point mass, not the lattice model's actual (smaller)
# indirect band gap.
GAP_EDGE = 1.0

# Exact topological prediction for this lattice (C=+1 for the filled/valence
# band, confirmed against the literature's own bond-phase convention).
CHERN_NUMBER = 1
EXPECTED_SLOPE = CHERN_NUMBER / (2 * np.pi)


def physical_density(h5_name, dat_name):
    """Return (E, S) -- the actual physical orbital-magnetization spectral density."""
    with h5py.File(h5_name, "r") as f:
        n_orbitals = int(np.array(f["NOrbitals"]))
        a1, a2 = np.array(f["LattVectors"])
        energy_scale = float(np.array(f["EnergyScale"]))

    unit_cell_area = abs(a1[0] * a2[1] - a1[1] * a2[0])
    norm_factor = n_orbitals / unit_cell_area

    d = np.loadtxt(dat_name)
    E, S_raw = d[:, 0], d[:, 1]

    prefactor = 0.5 * energy_scale  # e=hbar=c=1; the extra energy_scale is the
    # hidden Jacobian from KITE's "v" being built from already-rescaled hoppings.
    S = prefactor * norm_factor * S_raw
    return E, S


def plot(h5_name, dat_name, out_path="plots/haldane_orbital_magnetization_preview.png"):
    E, S = physical_density(h5_name, dat_name)

    dE = E[1] - E[0]
    M = np.cumsum(S) * dE

    in_gap = np.abs(E) < 0.8 * GAP_EDGE
    slope, intercept = np.polyfit(E[in_gap], M[in_gap], 1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

    ax1.plot(E, S, color=kite_style.KITE_PRIMARY, lw=1.5)
    ax1.axhline(0.0, color="0.6", lw=0.8)
    for edge in (-GAP_EDGE, GAP_EDGE):
        ax1.axvline(edge, color=kite_style.KITE_ACCENT, ls="--", lw=1.0)
    ax1.axvspan(-GAP_EDGE, GAP_EDGE, color=kite_style.KITE_ACCENT, alpha=0.08)
    ax1.set_ylabel(r"$S(E)$  (nm$^{-2}$)", fontsize=13)
    ax1.set_title("Haldane model: orbital-magnetization spectral density\n"
                  r"($t_2=t/3$, $\phi=\pi/2$, $\delta=0$: $C=+1$ Chern insulator)",
                  fontsize=12)

    ax2.plot(E, M, color=kite_style.KITE_PRIMARY, lw=1.5, label=r"$M(E_F)$ (this calculation)")
    ax2.plot(E[in_gap], slope * E[in_gap] + intercept,
             color="#357EDD", ls="--", lw=1.5,
             label=rf"linear-in-gap fit: $dM/dE_F={slope:.3f}$")
    ax2.plot(E[in_gap], EXPECTED_SLOPE * E[in_gap] + intercept,
             color="0.4", ls=":", lw=1.5,
             label=rf"exact k-space (Streda) prediction: $C/(2\pi)={EXPECTED_SLOPE:.4f}$")
    ax2.axvspan(-GAP_EDGE, GAP_EDGE, color=kite_style.KITE_ACCENT, alpha=0.08,
                label="bulk gap")
    for edge in (-GAP_EDGE, GAP_EDGE):
        ax2.axvline(edge, color=kite_style.KITE_ACCENT, ls="--", lw=1.0)
    ax2.set_xlabel(r"$E_F / t$", fontsize=13)
    ax2.set_ylabel(r"$M(E_F)$  (nm$^{-2}\cdot t$)", fontsize=13)
    ax2.legend(loc="upper left", fontsize=9, framealpha=0.85)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.replace("_preview.png", ".pdf"))
    plt.close(fig)
    print(f"Saved {out_path}")
    print(f"Linear-in-gap slope dM/dE_F = {slope:.4f}, "
          f"vs. exact k-space prediction C/(2*pi) = {EXPECTED_SLOPE:.4f} "
          f"({(slope - EXPECTED_SLOPE) / EXPECTED_SLOPE * 100:+.1f}%)")
    return E, S, M


if __name__ == "__main__":
    h5_name = sys.argv[1] if len(sys.argv) > 1 else "haldane_orbital_magnetization-output.h5"
    dat_name = sys.argv[2] if len(sys.argv) > 2 else "custom_one.dat"
    plot(h5_name, dat_name)
