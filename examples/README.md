!!! Info

    This file describes the contents of the [`#!bash kite/examples/`][examples-folder-github]-folder.

## KITE Scripts
This folder contains various scripts with examples showing the functionalities of the KITE library. To execute a script, run

``` bash
python3 example.py
../build/KITEx example-output.h5
../build/KITE-tools example-output.h5
```

These Python scripts return a HDF5-file. These files can then be passed to the KITEx-executable to perform the core calculation.
Finally, the raw output data can be analyzed using the KITE-tools-executable.
Depending on the type of calculation, various output files will be created in the folder where the code was executed from.

!!! Info
    A few flagship examples get a full step-by-step write-up elsewhere in these docs, under
    **In Depth Examples**: [Graphene][graphene-example], [Haldane][haldane-example],
    [Phosphorene][phosphorene-example], [Optical Conductivity][optical-conductivity-example], and
    [Spectral Function (ARPES)][spectral-function-example]. They are only listed briefly below, with a
    link to the deeper write-up, to avoid duplicating that material.

## `examples/*.py`

These are pure-KITE scripts (using `kite.lattice.Lattice`, `src/kite/lattice.py`) — no pybinding dependency required.

**Density of states (DOS)**

| Script | What it demonstrates |
| --- | --- |
| `dos_graphene.py` | DOS of monolayer graphene. See the [in-depth Graphene write-up][graphene-example]. |
| `dos_checkerboard_lattice.py` | DOS of a 2D checkerboard lattice. |
| `dos_cube.py` | DOS of a simple cubic lattice (3D). |
| `dos_cubic_lattice_twisted_bc.py` | DOS of a simple cubic lattice with random twisted boundary conditions. |
| `dos_square_lattice.py` | DOS of a square lattice, automatic energy rescaling. |
| `dos_square_lattice_disorder.py` | DOS of a square lattice with on-site (`kite.Disorder`) disorder. |
| `dos_square_lattice_twisted_bc.py` | DOS of a square lattice with twisted boundary conditions. |
| `dos_t_symmetric_cubic_weyl_sm.py` | DOS of a time-reversal-symmetric cubic Weyl semimetal (two orbitals per site) with Anderson disorder [Pixley, Goswami & Das Sarma, PRB 93, 085103 (2016)]. |
| `dos_fu_kane_mele_model.py` | DOS of the Fu-Kane-Mele topological insulator on the diamond lattice [Fu, Kane & Mele, PRL 98, 106803 (2007)]. |
| `dos_vacancies.py` | Honeycomb lattice with structural vacancy disorder (`kite.StructuralDisorder`) at two different concentrations, one per sublattice. |
| `dos_dccond_square_lattice.py` | DOS and DC conductivity of a square lattice, combined in one script. |
| `dos_dccond_haldane.py` | DOS and DC (xy) conductivity of the Haldane model, with uniform on-site disorder on both sublattices. |

**Optical / DC conductivity**

| Script | What it demonstrates |
| --- | --- |
| `hbn_optcond2_vacancies.py` | Optical conductivity of hexagonal boron nitride with vacancy disorder. |
| `optcond_t_symmetric_cubic_weyl_sm.py` | Linear optical conductivity (\(\sigma_{xx}\)) of the T-symmetric cubic Weyl semimetal; a heavily-parallelized example. See also the [Optical Conductivity write-up][optical-conductivity-example]. |
| `dccond_phosphorene.py` | Single-shot DC conductivity (`xx`/`yy`) of bilayer phosphorene. See the [in-depth Phosphorene write-up][phosphorene-example]. |
| `weyl_lt.py` | The same Weyl-semimetal optical-conductivity calculation as `optcond_t_symmetric_cubic_weyl_sm.py`, built using KITE's own `kite.lattice.Lattice` API rather than pybinding — compare directly against `pybinding/weyl_pb.py`. |

!!! warning "Known issue"
    `dccond_phosphorene.py` currently has a pre-existing, unrelated bug: line 145 references an undefined
    name `npoints` (the variable actually defined a few lines above is `num_points`), so running the script
    as-is raises a `NameError`. This is independent of anything else described on this page.

**Local density of states (LDOS) and spectral function (ARPES)**

| Script | What it demonstrates |
| --- | --- |
| `ldos_graphene.py` | Local density of states of graphene at a set of chosen positions/orbitals. |
| `arpes_bilayer.py` | One-particle spectral function (ARPES-like) of bilayer graphene with Rashba spin-orbit coupling, with uniform on-site disorder applied identically across all 8 sublattices. See [Disorder: adding the same disorder to several sublattices at once][disorder-list-shortcut] for how this script's disorder setup can be simplified. |
| `arpes_cubic.py` | One-particle spectral function of a simple cubic lattice with on-site disorder. |
| `arpes_tmd.py` | One-particle spectral function of a monolayer transition-metal dichalcogenide (3-band `dz2`/`dxy`/`dx2-y2` model), built from `kite.repository.group6_tmd.monolayer_3band()`. See the [in-depth Spectral Function write-up][spectral-function-example]. |
| `altermagnet_arpes.py` | Spin-resolved ARPES on a minimal 2D d-wave altermagnet (square lattice, Néel exchange + sublattice-swapped anisotropic NNN hopping): spin-up and spin-down Fermi pockets come out exactly 90°-rotated copies of each other, the defining altermagnet signature (zero net moment, non-relativistic spin splitting). Post-process with `process_altermagnet_arpes.py`. See the [in-depth write-up][altermagnet-example], including the full Bloch Hamiltonian. |
| `piflux_ldos_map.py` | Real-space local-density-of-states map (`calculation.ldos_map`) around a single vacancy in a π-flux (Dirac) square lattice, at the Dirac-point energy. Post-process with `process_piflux_ldos.py`. See the [in-depth write-up][markov-maps-example]. |
| `weyl_spectral_map.py` | Momentum-space spectral-function map (`calculation.spectral_map`) of a 3D Weyl semimetal, showing the Weyl-node iso-energy ring contours. Post-process with `process_weyl_spectral.py`. See the [in-depth write-up][markov-maps-example]. |

**Custom two-operator (Vertex) correlation functions**

| Script | What it demonstrates |
| --- | --- |
| `shinada_single.py` | A 3-orbital (`d`, `px`, `py`) "cuprate-like" lattice with complex Rashba-like hoppings, evaluated with KITE's `custom.Vertex`/`calculation.custom_singleshot_two` machinery at a fixed list of energies. |
| `shinada_fs.py` | The same lattice and Vertex construction as `shinada_single.py`, but swept over an energy grid at finite temperature via `calculation.custom_two`. |
| `kane_mele_spin_hall.py` | The Kane-Mele model's quantized spin Hall conductivity, computed via `custom.Vertex` + `calculation.custom_two` (spin current vertex `A = (1/2){v_x, s_z}`, `s_z` built from `add_orbital_index`/`add_orbital_coupling`). Shows a genuinely flat plateau at σ ≈ −2.02 across the bulk gap. Post-process with `kane_mele_spin_hall_process.py`. See the [in-depth write-up][custom-vertex-example]. |
| `kane_mele_spin_hall_disorder.py` | Extension of `kane_mele_spin_hall.py`: sweeps onsite Anderson disorder strength and shows the spin Hall plateau eroding once disorder approaches/exceeds the bulk gap scale. See the [in-depth write-up][custom-vertex-example]. |

**Custom local potential**

| Directory | What it demonstrates |
| --- | --- |
| `custom_local_potential/` | A worked example of KITE's custom local (site- and orbital-dependent) potential mechanism. See the [dedicated write-up][custom-local-potential-example] for the full explanation. |

**Utility / plotting helper scripts**

The folder also contains a handful of non-calculation helper scripts that are not examples in their own
right: `process_cond.py`, `process_conductivities.py`, `process_statistics.py`, and `cond_sum.py`
(Chebyshev-moment summation/post-processing helpers used by some of the scripts above), plus
`myaux.py`, `kite_style.py`, and `kite.mplstyle` (shared plotting/style helpers).

## `examples/pybinding/*.py`

The scripts in this subfolder still depend on [Pybinding](https://docs.pybinding.site/), which is now an
optional extra for KITE (`pip install -e ".[pybinding]"`), not a core dependency. See
[`examples/pybinding/README.md`][pybinding-readme-github] for the up-to-date, detailed rationale; in short:

| Script | Why it still needs pybinding |
| --- | --- |
| `dos_mixed_disorder.py`, `dos_on_site_disorder.py`, `dos_optcond_gaussian_disorder.py` | Use `pybinding.repository.graphene`, a built-in lattice template, instead of defining the lattice by hand. |
| `dos_twisted_bilayer.py` | Uses pybinding's own `Model`/`kpm`/`translational_symmetry` workflow directly, not just lattice construction — a fundamentally different workflow from KITE's. |
| `weyl_pb.py` | Builds its lattice with `pb.Lattice(...)` for direct comparison against the equivalent top-level `weyl_lt.py`, which uses KITE's own lattice API instead. |

`run_all_examples.py` in this folder runs and cleans up all of the pybinding-dependent examples above.

## `examples/paper/`

This folder reproduces the figures from KITE's own paper (`rsos.191809.pdf`), organized into one
subdirectory per paper section. A one-line summary of each is enough here; see
[Examples from KITE's Paper][paper-examples-page] for more, and browse each subdirectory's `Input`/`Output`
folders directly for the full scripts.

| Subdirectory | Summary |
| --- | --- |
| `Section_4_A_twisted_bilayer` | DOS of twisted bilayer graphene, relaxed and unrelaxed, at several energy resolutions. |
| `Section_4_B_topology` | Quantum anomalous Hall effect: optical and DC conductivity of a Kane-Mele honeycomb lattice with Rashba SOC, exchange terms, and Gaussian disorder. |
| `Section_4_D_electronic_structure_I` | ARPES (one-particle spectral function) of bilayer graphene with Rashba SOC and vacancy disorder. |
| `Section_4_D_electronic_structure_II` | Local density of states (LDOS) of the Haldane model with vacancy disorder. |
| `Section_4_E_spintronics` | Time-evolution of a Gaussian wave-packet in a Kane-Mele/Rashba honeycomb lattice, with and without magnetic (Anderson) disorder. |
| `Section_4_F_Magnetic_Field` | DOS of bilayer phosphorene in a magnetic field, at large system sizes. |
| `Section_4_G_molecules` | DOS of individual polycyclic aromatic hydrocarbon molecules (benzene through pentacene/oligoacenes) with on-site disorder. |
| `Section_6_benchmarks` | Performance benchmarking scripts (DOS of graphene with nearest- and next-nearest-neighbor hoppings) used to generate the paper's scaling/benchmark plots. |

## Running everything at once

All the results can be generated automatically by running

``` bash
python3 run_all_examples.py
```

After running this command, all the examples will be executed. This can take several minutes.
Besides the output files, like *name-dos.dat*, plots will be given for the *DOS*, *optical conductivity* and *DC conductivity*.

To clean up the folder after running all the examples, execute the following commands

``` bash
python3
>>> import run_all_examples
>>> run_all_examples.clean()
>>> exit()
```

[examples-folder-github]: https://github.com/quantum-kite/kite-v2/tree/master/examples
[pybinding-readme-github]: https://github.com/quantum-kite/kite-v2/tree/master/examples/pybinding
[paper-examples-page]: paper.md
[custom-local-potential-example]: custom_local_potential.md
[disorder-list-shortcut]: ../disorder.md#adding-the-same-disorder-to-several-sublattices-at-once
[markov-maps-example]: ../examples/markov_local_maps.md
[custom-vertex-example]: ../examples/custom_vertex_operators.md
[altermagnet-example]: ../examples/altermagnet_arpes.md
[graphene-example]: ../examples/graphene.md
[haldane-example]: ../examples/haldane.md
[phosphorene-example]: ../examples/phosphorene.md
[optical-conductivity-example]: ../examples/optical_conductivity.md
[spectral-function-example]: ../examples/spectral_function.md
