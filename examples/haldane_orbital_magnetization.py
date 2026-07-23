""" Haldane model: real-space spectral approach to orbital magnetization

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    This example demonstrates KITE's rank-ONE custom-vertex machinery,
    Calculation.custom_one(), on the paradigmatic single-spin Chern insulator:
    the Haldane model. Where custom_two() builds a DOUBLE-Chebyshev correlator
    Tr[T_m(H) B T_n(H) A] of two operators (the Kubo-Bastin trace behind the
    conductivity / spin-Hall responses, see kane_mele_spin_hall.py and the
    custom_vertex_operators.md doc page), custom_one() builds a SINGLE-Chebyshev
    moment vector

        mu_n = Tr[ T_n(H) A ]                                     (n = 0..M-1)

    of ONE operator A. Reconstructed with the standard KPM kernel, mu_n gives
    the operator-weighted spectral density  A(E) = Tr[ A delta(E - H) ], whose
    running integral up to a Fermi level mu returns the ground-state expectation
    Tr[ A theta(mu - H) ] = sum over occupied states of <k|A|k>.

    Here A is the (symmetrized) orbital magnetic moment operator, so A(E) is an
    energy-resolved orbital-magnetization spectral function and its integral up
    to E_F is the total orbital magnetization -- exactly the "real-space spectral
    approach to orbital magnetization" of Vidarte et al., Phys. Rev. B,
    doi 10.1103/bhg8-c3rw (arXiv:2512.01575).

    The orbital-magnetization vertex
    --------------------------------
    Operator strings apply RIGHT-TO-LEFT: "rx.vy" = r_x v_y (v_y acts first).
    The vertex used below is the RAW, un-symmetrized commutator form

        A = rx.vy - ry.vx = x H y - y H x

    with plain real coefficients (+1, -1) -- no manually-inserted factor of i.

    KITE's raw "vx"/"vy" building block is [H, r] (no i -- see custom_one's own
    docstring in src/kite/__init__.py for the derivation from the C++ source),
    so A above is genuinely ANTI-Hermitian (it's built from an odd number -- one
    -- of "v" tokens per term). Tr[T_n(H) A] is therefore purely IMAGINARY, not
    real: this is the correct, physical signal for this vertex, not an error.
    KITE auto-detects the "v"-token count from the vertex string itself
    (Calculation.custom_one counts them and writes NumVelocities to the HDF5;
    Src/Simulation/Custom/SimulationRankOne.cpp's store_custom_one then keeps
    the imaginary part instead of the real part when that count is odd -- the
    same parity bookkeeping Tools/Gamma2D.cpp already used for the built-in
    conductivity_dc/conductivity_optical family). No manual insertion of i is
    needed or wanted: Chebyshev moments don't care whether the traced operator
    is Hermitian or anti-Hermitian, so the correction belongs in the
    Hermitization/reconstruction step, not baked into the vertex itself.

    Relating A to the actual orbital magnetization operator, and the sign of the
    physical prefactor that belongs in POST-PROCESSING (not here): the reference
    (Vidarte et al., Eq. 1-2) defines
        M_hat_z = -(i e)/(2 hbar c Area) * (x H y - y H x) = -(e/(2 hbar c Area)) * A,
    i.e. a NEGATIVE real prefactor multiplying this exact vertex (see
    process_haldane_orbital_magnetization.py for where that prefactor and the
    Chebyshev-rescaling Jacobian 1/energy_scale are applied).

    The Haldane model and why this phase is topological
    ---------------------------------------------------
    Honeycomb lattice, two sublattices A/B. Nearest-neighbour hopping -t is
    spin-diagonal charge hopping; the next-nearest-neighbour hopping carries the
    complex Haldane phase, +phi on the A-A bonds and -phi on the B-B bonds. We
    use t2 = t/3 and the canonical flux phi = pi/2, kept EXPLICIT here because it
    is a genuine tunable of the Haldane construction (unlike kane_mele_spin_hall.py,
    which hardcodes the equivalent +/- i on the NNN bonds, i.e. phi = pi/2 too).

    The two Dirac points acquire opposite-sign Haldane masses of magnitude
    3*sqrt(3)*t2*sin(phi); a sublattice on-site asymmetry delta (onsite=(delta,-delta))
    adds a valley-independent Semenoff mass. The gap closes at
    |delta| = 3*sqrt(3)*t2*sin(phi), and the phase is topological (Chern number
    C = +/-1) for |delta| below that. Here 3*sqrt(3)*t2*sin(phi) = sqrt(3) ~ 1.73,
    so delta = 0 (pure Haldane) sits deep in the C = +/-1 phase -- a clean, gapped
    Chern insulator. We therefore keep delta = 0.

    We also DROP the small on-site disorder of the original script for this first
    clean demonstration (matching the precedent set by kane_mele_spin_hall.py,
    which keeps its main example clean and puts disorder in a separate extension).

    Boundaries: OPEN. The position operator r_x, r_y entering A is only
    well-defined with open boundaries (r is unbounded under periodic wrapping),
    so open boundaries are mandatory for this position-operator calculation --
    this also matches the original reference script.

    Units: energy in units of hopping |t| = 1, length in nm.
    Last updated: 22/07/2026
"""

__all__ = ["haldane", "main"]

import kite
import numpy as np
from kite import custom
from kite import lattice as latt


def haldane(t=1.0, phi=0.5 * np.pi, delta=0.0):
    """Return the Haldane honeycomb lattice (2 sublattices A, B).

    t2 = t/3, complex NNN phase +phi on A-A and -phi on B-B, on-site
    asymmetry (delta, -delta). delta = 0 => pure Haldane Chern insulator.
    """
    a = 0.24595   # [nm] unit-cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    t2 = t / 3.0
    phase = np.exp(1j * phi)

    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, 0.5 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)

    # Order fixes alias_id: A = 0, B = 1
    lat.add_sublattices(
        ('A', [0, -a_cc / 2],  delta),
        ('B', [0,  a_cc / 2], -delta),
    )

    lat.add_hoppings(
        # Nearest-neighbour charge hopping -t
        ([0, 0],  'A', 'B', -t),
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t),
        # Next-nearest-neighbour Haldane hopping: +phi on A-A, -phi on B-B
        ([1, 0],  'A', 'A', t2 * phase),
        ([0, -1], 'A', 'A', t2 * phase),
        ([-1, 1], 'A', 'A', t2 * phase),
        ([1, 0],  'B', 'B', t2 * np.conj(phase)),
        ([0, -1], 'B', 'B', t2 * np.conj(phase)),
        ([-1, 1], 'B', 'B', t2 * np.conj(phase)),
    )
    return lat


def main(t=1.0, phi=0.5 * np.pi, delta=0.0, moments=512, num_random=40,
         output_file="haldane_orbital_magnetization-output.h5"):
    """Prepare the input file for KITEx (orbital magnetization via custom_one)."""
    lattice = haldane(t, phi, delta)

    nx = ny = 2
    lx = ly = 64

    # OPEN boundaries: the position operator in the vertex requires them.
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=['open', 'open'],
        is_complex=True,
        precision=1,
        spectrum_range=[-4.0, 4.0],
    )

    calculation = kite.Calculation(configuration)

    # A = x*H*y - y*H*x, the raw commutator form -- no manually-inserted i (see
    # module docstring). One velocity ("v") token per term -> odd count -> KITE
    # auto-detects this and keeps the imaginary part in store_custom_one.
    A = custom.Vertex(moments, [[1.0, "rx.vy"],
                                 [-1.0, "ry.vx"]])

    calculation.custom_one(
        stream_=A,
        num_random_=num_random,
        num_disorder_=1,
    )

    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # Run:   ../build/KITEx haldane_orbital_magnetization-output.h5
    # There is NO KITE-tools step for custom_one -- the Chebyshev moment vector
    # /Calculation/CustomOne/Gamma is reconstructed into an energy-resolved
    # orbital-magnetization spectral density in Python; see
    # process_haldane_orbital_magnetization.py.
    return output_file


if __name__ == "__main__":
    main()
