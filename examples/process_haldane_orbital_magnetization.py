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

    Physical normalization
    ----------------------
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
       are built directly from the hopping values as exported to the HDF5 --
       which src/kite/__init__.py divides by EnergyScale before KITEx runs.
       So the vertex KITE actually evaluates is
       A_computed = A_physical / EnergyScale, and rho_A(E) (physical) needs an
       extra factor of EnergyScale to compensate.

    With e = hbar = c = 1: S(E) = 0.5 * EnergyScale * (NOrbitals/unit_cell_area) * S_raw(E).

    Full-range verification against the k-space modern theory
    -----------------------------------------------------------
    This real-space M(E_F) is checked against the k-space modern theory of
    orbital magnetization (Xiao, Shi & Niu; Thonhauser, Ceresoli, Vanderbilt &
    Resta) evaluated across the WHOLE plotted energy range, not just the bulk
    gap: for each band n(k),

        m_n(k)     = Im[ sum_{m!=n} Vx_nm Vy_mn ] / (E_n - E_m)
        Omega_n(k) = 2 Im[ sum_{m!=n} Vx_nm Vy_mn ] / (E_n - E_m)^2
        M(E_F) = sum_n integral d2k/(2*pi)^2 Theta(E_F - E_n(k))
                         * [ m_n(k) + (E_F - E_n(k)) * Omega_n(k) ]

    with Vx = U^dagger dH/dkx U, Vy = U^dagger dH/dky U in the band eigenbasis.
    kspace_full_range() below builds H(k) and its analytic k-derivatives
    directly and vectorizes over the whole Brillouin-zone grid (batched
    numpy.linalg.eigh over the leading grid axes), then turns the per-k-point
    (E_n, m_n, Omega_n) into M(E_F) at arbitrary E_F via a sort + cumulative
    sum (the E_F integral is a step function of E_n, so this avoids
    re-summing over the whole BZ for every E_F).
"""

import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt

import kite_style
from haldane_orbital_magnetization import haldane

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


def _hopping_terms_and_onsite(lattice):
    """Return (terms, onsite): one-directional hopping terms as
    (t, to_id, from_id, total_vector), and the (n, n) onsite matrix.

    total_vector = R_cartesian + (d_to - d_from), matching the exact phase
    convention of kite.visualize.hamiltonian_k -- see that function's
    docstring for the derivation (minus sign in the exponent).
    """
    n = lattice.nsub
    onsite = np.zeros((n, n), dtype=complex)
    for sub in lattice.sublattices.values():
        onsite[sub.alias_id, sub.alias_id] = np.asarray(sub.energy)[0, 0]

    vectors_matrix = np.array(lattice.vectors, dtype=float)
    terms = []
    for hop in lattice.hoppings.values():
        t = np.asarray(hop.energy)[0, 0]
        d_a = np.asarray(lattice.sublattices[hop.from_sub].position, dtype=float)[:2]
        d_b = np.asarray(lattice.sublattices[hop.to_sub].position, dtype=float)[:2]
        for term in hop.terms:
            R = np.asarray(term.relative_index, dtype=float) @ vectors_matrix
            terms.append((t, term.to_id, term.from_id, R + (d_b - d_a)))
    return terms, onsite


def _hamiltonian_k_batch(lattice, KX, KY):
    """H(k), dH/dkx, dH/dky over a whole k-grid at once (shape (..., n, n)).

    Same gauge/sign convention as kite.visualize.hamiltonian_k, but built with
    analytic k-derivatives (phase = exp(-i k.total), so d(phase)/dk is exact,
    no finite differencing) and batched over the grid's leading axes via
    numpy broadcasting/matmul/eigh -- evaluating hamiltonian_k point-by-point
    in a Python loop does not scale to a Brillouin-zone grid fine enough for
    the integral below.
    """
    terms, onsite = _hopping_terms_and_onsite(lattice)
    n = lattice.nsub
    shape = KX.shape
    H0 = np.zeros(shape + (n, n), dtype=complex)
    dH0dx = np.zeros_like(H0)
    dH0dy = np.zeros_like(H0)
    for t, to_id, from_id, total in terms:
        phase = np.exp(-1j * (KX * total[0] + KY * total[1]))
        H0[..., to_id, from_id] += t * phase
        dH0dx[..., to_id, from_id] += t * (-1j * total[0]) * phase
        dH0dy[..., to_id, from_id] += t * (-1j * total[1]) * phase
    H = H0 + np.conj(np.swapaxes(H0, -1, -2)) + onsite
    dHdx = dH0dx + np.conj(np.swapaxes(dH0dx, -1, -2))
    dHdy = dH0dy + np.conj(np.swapaxes(dH0dy, -1, -2))
    return H, dHdx, dHdy


def kspace_full_range(lattice, Efs, Nk=80):
    """Full-range k-space M(E_F) via the modern theory (Xiao-Shi-Niu), vectorized.

    Nk=80 is already converged for this lattice (checked up to Nk=400: the
    in-gap slope and M(E_F) values away from band edges do not change beyond
    the 4th significant figure), since the gap stays wide (>=1) everywhere in
    the Brillouin zone -- no near-degenerate k-points that would need a finer
    grid to resolve.

    Returns an array the same shape as Efs.
    """
    b1, b2 = lattice.reciprocal_vectors()
    f1 = np.linspace(-0.5, 0.5, Nk, endpoint=False)
    f2 = np.linspace(-0.5, 0.5, Nk, endpoint=False)
    F1, F2 = np.meshgrid(f1, f2, indexing="ij")
    KX = F1 * b1[0] + F2 * b2[0]
    KY = F1 * b1[1] + F2 * b2[1]

    H, dHdx, dHdy = _hamiltonian_k_batch(lattice, KX, KY)
    E, U = np.linalg.eigh(H)
    Ud = np.conj(np.swapaxes(U, -1, -2))
    Vx = Ud @ dHdx @ U
    Vy = Ud @ dHdy @ U

    n_bands = E.shape[-1]
    BZ_area = abs(b1[0] * b2[1] - b1[1] * b2[0])
    weight = BZ_area / (Nk * Nk) / (2 * np.pi) ** 2

    E_all, C1_all, C2_all = [], [], []
    for n in range(n_bands):
        s_m = np.zeros(KX.shape)
        s_o = np.zeros(KX.shape)
        for m in range(n_bands):
            if m == n:
                continue
            num = Vx[..., n, m] * Vy[..., m, n]
            dE_nm = E[..., n] - E[..., m]
            # Sign fixed against the exact Streda-relation check (dM/dE_F =
            # C/(2*pi) in the gap): the textbook +Im[...] convention here
            # matches this vertex/electron-charge sign; -Im[...] gave the
            # exact negative of the expected slope.
            s_m += np.imag(num) / dE_nm
            s_o += 2 * np.imag(num) / dE_nm ** 2
        En = E[..., n]
        E_all.append(En.ravel())
        C1_all.append((s_m - En * s_o).ravel())  # coefficient of E_F^0
        C2_all.append(s_o.ravel())                # coefficient of E_F^1

    E_all = np.concatenate(E_all)
    order = np.argsort(E_all)
    E_sorted = E_all[order]
    cum1 = np.concatenate([[0.0], np.cumsum(np.concatenate(C1_all)[order]) * weight])
    cum2 = np.concatenate([[0.0], np.cumsum(np.concatenate(C2_all)[order]) * weight])

    idx = np.searchsorted(E_sorted, Efs, side="right")
    return Efs * cum2[idx] + cum1[idx]


def plot(h5_name, dat_name, out_path="plots/haldane_orbital_magnetization_preview.png"):
    E, S = physical_density(h5_name, dat_name)

    dE = E[1] - E[0]
    M = np.cumsum(S) * dE

    in_gap = np.abs(E) < 0.8 * GAP_EDGE
    slope, intercept = np.polyfit(E[in_gap], M[in_gap], 1)

    lattice = haldane()
    M_kspace = kspace_full_range(lattice, E)

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

    ax2.plot(E, M_kspace, color="0.4", ls=":", lw=2.0,
             label="k-space modern theory (full range)")
    ax2.plot(E, M, color=kite_style.KITE_PRIMARY, lw=1.5, label=r"$M(E_F)$ (real-space, this calculation)")
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
    rmse = np.sqrt(np.mean((M - M_kspace) ** 2))
    scale = np.abs(M_kspace).max()
    print(f"Full-range real-space vs. k-space RMSE = {rmse:.4f} "
          f"({rmse / scale * 100:.1f}% of peak |M|)")
    return E, S, M, M_kspace


if __name__ == "__main__":
    h5_name = sys.argv[1] if len(sys.argv) > 1 else "haldane_orbital_magnetization-output.h5"
    dat_name = sys.argv[2] if len(sys.argv) > 2 else "custom_one.dat"
    plot(h5_name, dat_name)
