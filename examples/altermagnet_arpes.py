""" 2D d-wave altermagnet (minimal square-lattice model): spin-split ARPES pockets

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Minimal square-lattice model of a d-wave altermagnet: two magnetic
    sublattices A, B (a checkerboard arrangement, A at [0,0] and B at
    [0.5,0.5] within the unit cell), spin-full (4 sublattices: Aup, Bup,
    Adn, Bdn). Altermagnets have zero net magnetic moment (A and B carry
    opposite Neel exchange fields, and up/down spin see the opposite sign of
    that field) yet a NON-relativistic spin splitting of the bands, because
    the exchange field is combined with an anisotropic (here, NNN) hopping
    that swaps sign between A and B. The result is d-wave momentum-dependent
    spin splitting: eps_up(kx,ky) = eps_dn(ky,kx), i.e. the up- and
    down-spin Fermi pockets are 90-degree ROTATED copies of one another
    rather than being related by inversion/time-reversal as in a
    conventional antiferromagnet. See Smejkal, Sinova, Jungwirth,
    PRX 12, 040501 (2022).

    Model
    -----
    a1=[1,0], a2=[0,1] (square lattice, a=1), hopping t=1, Neel exchange
    J=0.8 (onsite, alternating sign on A/B and on up/down), anisotropic
    NNN delta=0.4:
      - NN A-B hopping -t (checkerboard bonds only, i.e. only between
        opposite sublattices -- see the geometry check in `main()`/module
        docstring below), identical for both spins (no SOC in this model).
      - NNN (same-sublattice) hopping is spin-independent but
        SUBLATTICE-SWAPPED: on A, -delta along x and +delta along y; on B,
        +delta along x and -delta along y. This sign swap between A and B
        is what makes the model d-wave (90-degree rotated) rather than the
        isotropic (s-wave, non-split) case.
    All hoppings are real; the up and down spin blocks are decoupled
    (block-diagonal 2x2 + 2x2), with splitting arising purely from
    combining the (spin-dependent) onsite Neel field with the
    (spin-independent, sublattice-dependent) anisotropic NNN hopping.

    Model verification (see docstring of `_verify_model()` and the printed
    output when this file is run directly): a direct 4x4 Bloch-Hamiltonian
    diagonalization (built from the exact same hopping list fed to KITE,
    not a textbook formula) confirms, on a k-grid:
      (a) zero net magnetic moment: Tr H(k) = 0 EXACTLY at every k (not
          just on average) -- both spin blocks are separately traceless.
      (b) eps_up(kx,ky) == eps_dn(ky,kx) to machine precision (the 90-degree
          altermagnet signature).
      (c) eps_up != eps_dn at generic k, but the two coincide exactly on
          the kx=ky and kx=-ky diagonals (nodal lines of the d-wave
          splitting) -- this is a real feature of the model, not a
          verification bug: on both diagonals, d(k) = 2*delta*(cos(ky) -
          cos(kx)) vanishes identically, so the two spin blocks become
          degenerate there.

    Calculation
    -----------
    Two independent KITE runs (KITE only allows a single `arpes()` request
    per Configuration/output file -- see `Calculation.arpes` in
    src/kite/__init__.py), one per spin, with orbital weights selecting
    only that spin's two sublattices ([1,1,0,0] for up, [0,0,1,1] for
    down). Both use the SAME 2D k-grid spanning the full square-lattice
    Brillouin zone [-pi,pi] x [-pi,pi] (a1=a2=1, so the reciprocal lattice
    vectors are exactly 2*pi*[1,0] and 2*pi*[0,1] and this range is the
    full BZ) -- a 2D grid rather than a 1D high-symmetry k-path was chosen
    specifically because the thing we want to show (a 90-degree rotation
    of the whole pocket shape) is far more visible as a 2D map than as a
    line cut; see `process_altermagnet_arpes.py` for the fixed-energy
    slice/reshape used to turn the ARPES A(k,E) output back into 2D maps.

    Units: energy in units of hopping |t| = 1, length in units of lattice
    constant |a| = 1.
    Last updated: 22/07/2026
"""

__all__ = ["altermagnet", "main"]

import kite
import numpy as np
from kite import lattice as latt


def altermagnet(t=1.0, J=0.8, delta=0.4):
    """Return the spin-full 2D d-wave altermagnet lattice (4 sublattices).

    Sublattice add order fixes alias_id: Aup=0, Bup=1, Adn=2, Bdn=3.
    """
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.0, 1.0])

    lat = latt.Lattice(a1=a1, a2=a2)

    # Order matters: alias_id 0,1,2,3. Onsite energy = Neel exchange field
    # (opposite sign on A vs B, opposite sign on up vs down spin).
    lat.add_sublattices(
        ('Aup', [0.0, 0.0], +J),
        ('Bup', [0.5, 0.5], -J),
        ('Adn', [0.0, 0.0], -J),
        ('Bdn', [0.5, 0.5], +J),
    )

    lat.add_hoppings(
        # --- NN A-B checkerboard bonds, amplitude -t, both spins ---
        # (KITE auto-adds the Hermitian conjugate for each of these, so
        # this list alone reproduces exactly the 4 A-B nearest-neighbour
        # bonds around every A site -- see `_verify_model()`'s geometric
        # distance check, which confirms these 4 bonds are equidistant).
        ([0, 0], 'Aup', 'Bup', -t),
        ([-1, 0], 'Aup', 'Bup', -t),
        ([0, -1], 'Aup', 'Bup', -t),
        ([-1, -1], 'Aup', 'Bup', -t),
        ([0, 0], 'Adn', 'Bdn', -t),
        ([-1, 0], 'Adn', 'Bdn', -t),
        ([0, -1], 'Adn', 'Bdn', -t),
        ([-1, -1], 'Adn', 'Bdn', -t),

        # --- Anisotropic NNN, spin-independent, sublattice-swapped ---
        # (each entry's Hermitian-conjugate/negative-R partner is added
        # automatically by KITE, exactly as for the Kane-Mele NNN terms in
        # examples/kane_mele_spin_hall.py).
        ([1, 0], 'Aup', 'Aup', -delta),
        ([0, 1], 'Aup', 'Aup', +delta),
        ([1, 0], 'Bup', 'Bup', +delta),
        ([0, 1], 'Bup', 'Bup', -delta),
        ([1, 0], 'Adn', 'Adn', -delta),
        ([0, 1], 'Adn', 'Adn', +delta),
        ([1, 0], 'Bdn', 'Bdn', +delta),
        ([0, 1], 'Bdn', 'Bdn', -delta),
    )
    return lat


def _verify_model(t=1.0, J=0.8, delta=0.4, n=61):
    """Direct 4x4 Bloch-Hamiltonian check of the exact hopping list above.

    Builds H(k) bond-by-bond (including the KITE-implicit Hermitian
    conjugate) rather than trusting a pre-derived textbook formula, and
    checks: (a) Tr H(k) == 0 at every k, (b) eps_up(kx,ky) == eps_dn(ky,kx)
    to machine precision, (c) eps_up != eps_dn at a generic k, but the two
    coincide on the kx=ky / kx=-ky diagonals. Raises AssertionError if any
    check fails -- this must pass before the model is trusted for a KITEx
    run.
    """
    # sublattice order: 0=Aup,1=Bup,2=Adn,3=Bdn
    onsite = np.array([+J, -J, -J, +J])

    nn_bonds = [
        ([0, 0], 0, 1, -t), ([-1, 0], 0, 1, -t),
        ([0, -1], 0, 1, -t), ([-1, -1], 0, 1, -t),
        ([0, 0], 2, 3, -t), ([-1, 0], 2, 3, -t),
        ([0, -1], 2, 3, -t), ([-1, -1], 2, 3, -t),
    ]
    nnn_bonds = [
        ([1, 0], 0, 0, -delta), ([0, 1], 0, 0, +delta),
        ([1, 0], 1, 1, +delta), ([0, 1], 1, 1, -delta),
        ([1, 0], 2, 2, -delta), ([0, 1], 2, 2, +delta),
        ([1, 0], 3, 3, +delta), ([0, 1], 3, 3, -delta),
    ]
    bonds = nn_bonds + nnn_bonds

    def bloch_H(kx, ky):
        H = np.diag(onsite).astype(complex)
        for R, i, j, val in bonds:
            phase = np.exp(1j * (kx * R[0] + ky * R[1]))
            H[i, j] += val * phase
            H[j, i] += np.conj(val * phase)  # KITE auto adds this bond
        return H

    # Geometry check: the 4 NN A-B bonds are equidistant (true nearest neighbours).
    A_pos = np.array([0.0, 0.0])
    B_local = np.array([0.5, 0.5])
    dists = []
    for R in [[0, 0], [-1, 0], [0, -1], [-1, -1]]:
        Bpos = R[0] * np.array([1.0, 0.0]) + R[1] * np.array([0.0, 1.0]) + B_local
        dists.append(np.linalg.norm(Bpos - A_pos))
    assert np.allclose(dists, dists[0]), f"A-B bond distances not equal: {dists}"

    kxs = np.linspace(-np.pi, np.pi, n)
    kys = np.linspace(-np.pi, np.pi, n)
    eps_up = np.zeros((n, n, 2))
    eps_dn = np.zeros((n, n, 2))
    trace_max = 0.0
    for i, kx in enumerate(kxs):
        for j, ky in enumerate(kys):
            H = bloch_H(kx, ky)
            assert np.allclose(H, H.conj().T), "H(k) not Hermitian"
            trace_max = max(trace_max, abs(np.trace(H).real))
            Hup = H[np.ix_([0, 1], [0, 1])]
            Hdn = H[np.ix_([2, 3], [2, 3])]
            eps_up[i, j] = np.sort(np.linalg.eigvalsh(Hup))
            eps_dn[i, j] = np.sort(np.linalg.eigvalsh(Hdn))

    assert trace_max < 1e-10, f"(a) net moment check failed: max|Tr H(k)|={trace_max:.3e}"

    rot_err = np.max(np.abs(eps_up - np.transpose(eps_dn, (1, 0, 2))))
    assert rot_err < 1e-8, f"(b) 90-deg rotation check failed: max err={rot_err:.3e}"

    # (c) generic k off both diagonals must split; on-diagonal k must NOT.
    kx_g, ky_g = 1.1, 0.3
    Hup_g = bloch_H(kx_g, ky_g)[np.ix_([0, 1], [0, 1])]
    Hdn_g = bloch_H(kx_g, ky_g)[np.ix_([2, 3], [2, 3])]
    split_generic = np.max(np.abs(np.sort(np.linalg.eigvalsh(Hup_g))
                                   - np.sort(np.linalg.eigvalsh(Hdn_g))))
    assert split_generic > 1e-3, "(c) expected nonzero splitting at generic k"

    kx_d = 0.8
    for ky_d in (kx_d, -kx_d):
        Hup_d = bloch_H(kx_d, ky_d)[np.ix_([0, 1], [0, 1])]
        Hdn_d = bloch_H(kx_d, ky_d)[np.ix_([2, 3], [2, 3])]
        split_diag = np.max(np.abs(np.sort(np.linalg.eigvalsh(Hup_d))
                                    - np.sort(np.linalg.eigvalsh(Hdn_d))))
        assert split_diag < 1e-10, f"(c) unexpected splitting on diagonal at ky={ky_d}"

    print("Model verification passed:")
    print(f"  A-B NN distances: {dists} (all equal)")
    print(f"  (a) max|Tr H(k)| over grid = {trace_max:.3e}")
    print(f"  (b) max|eps_up(kx,ky) - eps_dn(ky,kx)| over grid = {rot_err:.3e}")
    print(f"  (c) generic k=({kx_g},{ky_g}) splitting = {split_generic:.4f}"
          f", on-diagonal splitting = {split_diag:.3e}")


def k_grid(n=41, kmax=np.pi):
    """Flat list of [kx, ky] points on an n x n grid spanning the full BZ."""
    kxs = np.linspace(-kmax, kmax, n)
    kys = np.linspace(-kmax, kmax, n)
    return [np.array([kx, ky]) for kx in kxs for ky in kys], kxs, kys


def main(t=1.0, J=0.8, delta=0.4, moments=512, n_kgrid=21,
         spin="up", output_file=None):
    """Prepare the input file for KITEx (ARPES for a single spin channel).

    spin : 'up' or 'down' -- selects the orbital weight vector; KITE only
    allows one arpes() request per Configuration, so up and down spin need
    two separate output files (two calls to this function).

    n_kgrid=21 (441 k-points total) was chosen after timing a KITEx run
    directly: at num_moments=512 on this 128x128x[2,2] lattice, each
    k-point costs ~0.18s of KITEx wall time on this machine, so a much
    denser grid (a 41x41 grid was tried first) becomes multi-minute per
    spin; 21x21 keeps each spin's KITEx run to roughly a minute while still
    resolving the pocket shapes well enough to see the 90-degree rotation.
    """
    lattice = altermagnet(t, J, delta)

    nx = ny = 2
    lx = ly = 128

    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=['periodic', 'periodic'],
        is_complex=True,
        precision=1,
        spectrum_range=[-6, 6],
    )

    calculation = kite.Calculation(configuration)

    k_vector, kxs, kys = k_grid(n=n_kgrid)

    if spin == "up":
        weight = [1, 1, 0, 0]
        output_file = output_file or "altermagnet_arpes_up-output.h5"
    elif spin == "down":
        weight = [0, 0, 1, 1]
        output_file = output_file or "altermagnet_arpes_down-output.h5"
    else:
        raise ValueError("spin must be 'up' or 'down'")

    calculation.arpes(
        k_vector=k_vector,
        weight=weight,
        num_moments=moments,
        num_disorder=1)

    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # Run:   ../build/KITEx <output_file>
    #        ../build/KITE-tools <output_file> --ARPES -K jackson -S
    return output_file, kxs, kys


if __name__ == "__main__":
    _verify_model()
    main(spin="up")
    main(spin="down")
