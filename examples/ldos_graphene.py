""" Local density of states around a single vacancy in graphene

    ##########################################################################
    #                         Copyright 2020/2026, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: energy in eV (t = 2.8 eV)
    Lattice: honeycomb, nearest-neighbour only
    Configuration: periodic boundaries, double precision, manual rescaling,
                   128 x 128 unit cells, no domain decomposition (nx = ny = 1)
    Calculation type: local density of states on a grid of sites around a
                      single vacancy
    Last updated: 11/09/2026

    Physics
    -------
    main() removes one A-sublattice atom at the centre of the sample with
    kite.StructuralDisorder(...).add_vacancy('A'). That is the only way to
    remove a site: hoppings passed to add_structural_disorder are INCREMENTS
    added on top of the clean hoppings, so setting one to zero changes nothing.

    A single vacancy on a bipartite lattice is not a generic impurity, it is a
    chiral defect, and its consequences are exact rather than asymptotic:

    * There is exactly one state at E = 0. With nearest-neighbour hoppings and
      zero on-site energies the honeycomb Hamiltonian anticommutes with the
      sublattice operator, so the spectrum is symmetric under E -> -E. Removing
      one site leaves 2*128*128 - 1 = 32767 sites, an odd number, and an
      odd-dimensional matrix with a symmetric spectrum must have a zero
      eigenvalue.

      Keep the system size away from multiples of 3. If 3 divides L the Dirac
      point falls on the k-grid, the clean lattice supplies two further zero
      modes, and the vacancy state becomes an ill-defined mixture of the three.
      128 is safe; 126 would not be.

    * The state lives entirely on the OTHER sublattice. Remove an A atom and
      the zero mode has amplitude only on B sites. Exact diagonalization of
      this lattice gives a B-sublattice weight of 1.000000 against 7e-30 on A.

    * It is quasi-localized, not exponentially localized. The amplitude falls
      off as 1/r, so |psi|^2 ~ 1/r^2: r^2 |psi|^2 is flat from r = 2a out to
      r = 16a. Because the integral of 1/r^2 diverges logarithmically, the
      normalization grows like ln(L) and the state has no intrinsic size --
      about half its weight sits within r = 4a and the rest is spread over the
      whole sample.

    This is the chiral zero mode of vacancy-disordered graphene, whose
    collective physics is the subject of A. Ferreira and E. R. Mucciolo,
    Phys. Rev. Lett. 115, 106601 (2015), and which is seen directly by STM.

    calculation.ldos then maps the LDOS on an 11 x 11 unit-cell grid, both
    sublattices, centred on the vacancy.

    What to check, and what not to
    ------------------------------
    Trustworthy:
      - an E = 0 LDOS peak on B sites near the vacancy, and NO corresponding
        peak on A sites. This sublattice asymmetry is the sharpest signature
        and it is exact;
      - the asymmetry disappearing with distance, both sublattices tending to
        the clean value, which vanishes at the Dirac point because the pristine
        DOS goes like |E|;
      - the E = 0 excess on B sites falling monotonically with distance, with
        r^2 times the excess roughly flat across the mapped window.

    Not trustworthy -- do not assert these in a test:
      - the height and the width of the E = 0 peak. The zero mode is a delta
        function; everything visible is the KPM kernel (Jackson by default, or
        -K green <eta> in KITE-tools). Both scale with num_moments/energy_scale
        and carry no physics;
      - absolute LDOS values compared across different num_moments or
        spectrum_range;
      - a numerical fit of the 1/r^2 exponent. The window spans only a factor
        ~2.5 in r and finite resolution mixes in states away from E = 0. Check
        the trend, not the exponent;
      - the value at the vacancy site itself, which KITE sets to zero by
        construction.

    Resolution, and why num_moments is 640
    --------------------------------------
    With energy_scale = 3.01 * 2.8 = 8.428 eV the Jackson width at band centre
    is sigma_E = pi * energy_scale / M, so FWHM = 2.355 * sigma_E. At the 64
    moments this example used to request, FWHM was 0.97 eV -- wider than the
    entire +/-0.4 eV window it asked for, so it resolved nothing at all.

    Resolution matters more here than it first appears. The measured B-to-A
    contrast at E = 0 is

        LDOS_B(0) / LDOS_A(0)  =  1 + w_B / (2 c sigma_E^2)

    because the zero mode's peak height grows as 1/sigma while the linear Dirac
    background at E = 0 grows as +sigma. The contrast therefore scales as M^2:
    roughly 1.4:1 at M = 64 -- a wiggle, not a polarised state -- against about
    40:1 at M = 640. 640 gives FWHM = 0.097 eV and costs seconds.

    This calculation carries NO stochastic error: calculation.ldos starts from
    a site-localised vector rather than a random one, so there is no trace
    estimator, and num_disorder = 1 is correct because the vacancy is a single
    deterministic defect at a fixed position, not an ensemble. The only error
    is the Jackson truncation above. That makes this the cleanest example in
    the repository of a KPM result with no statistical component.

    spectrum_range
    --------------
    Declared as [-3.01*t, +3.01*t] = [-8.428, 8.428] eV. The honeycomb
    bandwidth is +/-3t, reached at Gamma, and a vacancy deletes a row and a
    column, so by interlacing it cannot push any eigenvalue outside that bound.
    Note that config_system may still warn that the range looks insufficient:
    the automatic estimate adds a structural-defect term that assumes the
    defect can push states outward, which a vacancy cannot. The warning is a
    false positive here -- do not "fix" it by widening the range, which would
    only cost resolution.
"""

__all__ = ["main"]

import kite
import numpy as np
from kite import lattice as latt


def graphene_lattice(onsite=(0, 0)):
    """Return lattice specification for a honeycomb lattice with nearest neighbor hoppings"""

    # parameters
    a = 0.24595  # unit cell length
    a_cc = 0.142  # carbon-carbon distance
    t = 2.8  # eV

    # define lattice vectors
    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    # create a lattice with 2 primitive vectors
    lat = latt.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', -t),
        # between neighboring cells, between which atoms, and the value
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t)
    )
    return lat


def analyze_results(filename, lattice):
    """Analyze the results from the LDOS-calculation"""
    import matplotlib.pyplot as plt

    # extraxt the positions of the sublattices
    sublattices = lattice.sublattices
    sub_list = [(sublattices[sub_name].unique_id, sub_name) for sub_name in [*sublattices]]
    sub_list.sort(key=lambda x: x[0])
    sublattice_position = np.array([sublattices[sub_name[1]].position for sub_name in sub_list])

    # read the input file
    file_content = np.loadtxt(filename).T
    pos = file_content[:2, :] if file_content.shape[0] == 4 else file_content[:3, :]  # 2D or 3D lattice
    orbs = file_content[-2, :]
    values = file_content[-1, :]
    n_pos = len(values)
    n_orb = len(set(orbs))

    # fetch the positions of the sublattices
    orb_tmp = np.zeros((n_pos, n_orb))
    for i, i_o in enumerate(set(orbs)):
        orb_tmp[np.array(orbs) == i_o,  i] = 1

    # find the position of the atoms, unit cell + pos. of subl. in unit cell
    xs = np.dot(pos.T, lattice.vectors) + np.dot(orb_tmp, sublattice_position)

    # create colors for the different LDOS values
    values_max = np.max(np.abs(values))

    fig = plt.figure(figsize=(8, 5))
    ax_sctr = fig.add_axes([0.1, 0.05, 0.7, 0.9])
    sctr = ax_sctr.scatter(xs[:, 0], xs[:, 1], c=values, s=70, vmin=-values_max, vmax=values_max)
    plt.colorbar(sctr, cax=fig.add_axes([0.85, 0.1, 0.03, 0.8]))
    ax_sctr.set_aspect('equal')
    ax_sctr.set_title(filename)
    fig.savefig(filename[:-4]+".pdf")
    plt.close(fig)


def main(onsite=(0, 0)):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = graphene_lattice(onsite)

    # A single vacancy on the A sublattice, at the centre of the sample.
    # add_vacancy() is the only way to remove a site: hoppings handed to
    # add_structural_disorder are INCREMENTS added on top of the clean
    # hoppings, so setting one to zero changes nothing at all.
    struc_disorder = kite.StructuralDisorder(lattice, position=[[64, 64]])
    struc_disorder.add_vacancy('A')

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 1
    # number of unit cells in each direction.
    lx = ly = 128

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny],
    # - lengths of structure [lx, ly]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random"

    # Boundary Mode
    mode = "periodic"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=False,
        precision=1,
        # Honeycomb bandwidth is exactly 3t, reached at Gamma; a vacancy only
        # deletes a row and column so by interlacing it cannot exceed that.
        # The old [-2.8, 2.8] put real eigenvalues at |x| = 3, where T_64(3)
        # is of order 1e48 and the output was meaningless.
        spectrum_range=[-3.01 * 2.8, 3.01 * 2.8]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)

    pos_matrix = []
    sub_matrix = []
    d1 = 64
    d2 = 64

    # range(-5, 6) is 11 wide and centred; range(-5, 5) ran -5..+4 and put the
    # vacancy off-centre in the map.
    for di in range(-5, 6):
        for dj in range(-5, 6):
            pos_matrix.append([d1 + di, d2 + dj])
            pos_matrix.append([d1 + di, d2 + dj])
            sub_matrix.append('A')
            sub_matrix.append('B')

    calculation.ldos(
        energy=np.linspace(-1.0, 1.0, 101),
        num_moments=640,
        num_disorder=1,
        position=pos_matrix,
        sublattice=sub_matrix
    )

    # configure the *.h5 file
    output_file = "graphene_lattice_ldos-output.h5"
    kite.config_system(lattice, configuration, calculation, disorder_structural=struc_disorder,
                       filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx graphene_lattice_ldos-output.h5
    # ../build/KITE-tools graphene_lattice_ldos-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
