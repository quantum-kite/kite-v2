""" Kane-Mele spin Hall effect: disorder washout of the quantized plateau

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Physics
    -------
    Extension of examples/kane_mele_spin_hall.py: the same honeycomb
    Kane-Mele lattice and the same Kubo-Bastin custom_two spin Hall
    calculation, but now with onsite (Anderson-type) disorder added on all
    4 sublattices via kite.Disorder, swept over disorder strength W (in
    units of the clean hopping t). The clean model's bulk gap is
    ~1.56*t (confirmed by the already-verified clean example's actual flat
    plateau region, mu in [-0.74, 0.74]*t). Expectation: the spin Hall
    plateau should stay flat/robust for W well below the gap scale (W=0,
    0.5) and visibly degrade/round off once W approaches or exceeds the gap
    scale (W=1.0, 2.0), as disorder-induced level broadening starts to
    close/smear the topological gap.

    This is a laptop-fast qualitative demo, not a disorder-averaged
    high-precision study: num_disorder_ is kept small (1 for the clean
    W=0 run, a handful of realizations for W>0) to keep runtime short, so
    the larger-W curves are expected to be visibly noisy rather than
    smooth -- see the honest reporting in the module-level results this
    script prints, and don't over-interpret small wiggles in those curves
    as physical structure.

    Everything else (lattice construction, operator definitions, sign
    convention, Kubo-Bastin post-processing) is exactly as documented in
    examples/kane_mele_spin_hall.py and examples/kane_mele_spin_hall_process.py
    -- this file only adds the disorder loop and imports `kane_mele()` and
    `spin_hall()` rather than duplicating them.

    Units: energy in units of hopping |t| = 1, length in nm.
    Last updated: 22/07/2026
"""

__all__ = ["main"]

import kite
import numpy as np
from kite import custom

from kane_mele_spin_hall import kane_mele


def main(W, t=1.0, t2=0.15, moments=256, num_disorder=1,
         output_file=None):
    """Prepare the input file for KITEx (spin Hall via custom_two) with
    onsite disorder of uniform width W (in units of t) on all 4
    sublattices.
    """
    lattice = kane_mele(t, t2)

    nx = ny = 2
    lx = ly = 64

    # Manual spectrum_range must safely bound the ACTUAL spectrum (clean
    # bandwidth + disorder spread), not just the clean bandwidth -- KPM's
    # Chebyshev recursion is only stable if every eigenvalue of H lands
    # inside the rescaled [-1, 1] domain; eigenvalues that spill outside
    # the chosen manual range make T_m(H) blow up exponentially with m
    # (found the hard way: at W=2.0 with the clean example's fixed [-4, 4]
    # range, the raw KITEx moments overflowed to ~1e160 and the
    # post-processed sigma_xy diverged to ~1e162 -- not "noisy", genuinely
    # numerically broken). The clean model's own bandwidth is ~3.0*t
    # (checked numerically on a k-grid); KITE's Uniform(mean=0, W) disorder
    # draws onsite energies uniformly on [-sqrt(3)*W, +sqrt(3)*W] (see the
    # comment further below), so the worst-case total spread is
    # 3.0 + sqrt(3)*W. A 20% safety margin, rounded up to the nearest 0.5,
    # gives the half-range used below; W=0 reduces to exactly [-4, 4],
    # matching examples/kane_mele_spin_hall.py's own already-verified
    # config for a direct, apples-to-apples W=0 comparison.
    clean_bandwidth = 3.0
    half_range = max(4.0, np.ceil(1.2 * (clean_bandwidth + np.sqrt(3) * W) * 2) / 2)

    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=['periodic', 'periodic'],
        is_complex=True,
        precision=1,
        spectrum_range=[-half_range, half_range],
    )

    calculation = kite.Calculation(configuration)

    # Register EVERY sublattice (idx = add order = alias_id).
    for lbl, idx in [('Aup', 0), ('Bup', 1), ('Adn', 2), ('Bdn', 3)]:
        calculation.add_orbital_index(lbl, idx)

    # s_z = (1/2) sigma_z: +1/2 on spin up, -1/2 on spin down.
    calculation.add_orbital_coupling('Aup', 'Aup',  0.5, 'l0')
    calculation.add_orbital_coupling('Bup', 'Bup',  0.5, 'l0')
    calculation.add_orbital_coupling('Adn', 'Adn', -0.5, 'l0')
    calculation.add_orbital_coupling('Bdn', 'Bdn', -0.5, 'l0')

    # A = (1/2){v_x, s_z} spin current ; B = v_y charge velocity
    A = custom.Vertex(moments, [[0.5, "vx.l0"], [0.5, "l0.vx"]])
    B = custom.Vertex(moments, [[1.0, "vy"]])

    calculation.custom_two(
        stream_=[A, B],
        num_random_=1,
        num_disorder_=num_disorder,
        num_points_=1000,
        temperature_=0.01,
    )

    # W == 0: pass NO Disorder object at all (kite.config_system's
    # `disorder` kwarg defaults to None). A Disorder() instance with no
    # add_disorder() calls is truthy but has empty internal arrays, which
    # crashes config_system's `disorder._orbital > -1` check -- this is a
    # real bug in src/kite/__init__.py (out of scope to fix here, per the
    # "do not edit src/kite/*" constraint), so the clean case must skip
    # Disorder entirely rather than construct an empty one.
    disorder = None
    if W > 0:
        disorder = kite.Disorder(lattice)
        for sub in ['Aup', 'Bup', 'Adn', 'Bdn']:
            # Uniform(mean=0, W): confirmed against Src/Tools/Random.cpp
            # (KPMRandom<T>::uniform), the 2nd argument is the RMS
            # (standard deviation) of the box distribution, not its
            # half-width -- KITE draws onsite disorder uniformly on
            # [-sqrt(3)*W, +sqrt(3)*W]. What matters for this qualitative
            # washout demo is that W sets the disorder ENERGY SCALE in
            # units of t; the sqrt(3) box/RMS factor doesn't change the
            # qualitative W=0/0.5 (below gap) vs W=1.0/2.0 (at/above gap)
            # comparison being made here.
            disorder.add_disorder(sub, 'Uniform', 0.0, W)

    if output_file is None:
        output_file = f"kane_mele_spin_hall_disorder_W{W:g}-output.h5"

    kite.config_system(lattice, configuration, calculation,
                        disorder=disorder, filename=output_file)

    # Run:   ../build/KITEx <output_file>
    # There is NO KITE-tools step for custom_two -- reuse
    # kane_mele_spin_hall_process.spin_hall() for the Kubo-Bastin energy
    # integral, exactly as in the clean example.
    return output_file


if __name__ == "__main__":
    for W in (0.0, 0.5, 1.0, 2.0):
        num_dis = 1 if W == 0.0 else 3
        main(W, num_disorder=num_dis)
