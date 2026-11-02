## Usage


### Default usage

Its default usage is very simple:

``` bash
./build/KITE-tools archive.h5
```

where archive.h5 is the HDF file that stores the output of KITE. If KITE-tools does not find this output, it will return an error. The output of KITE-tools is a set of .dat files, one for each of the requested quantities. KITE-tools may be executed without any additional parameters; all the unspecified parameters required for the calculation will be set to sensible default values. At the moment, KITE-tools is able to compute the following quantities:

  * Local density of states (LDOS)
  * Angle-resolved photoemission spectroscopy (experimental) (ARPES)
  * Density of states (DOS)
  * DC conductivity (CondDC)
  * Optical conductivity (CondOpt)
  * Second-order optical conductivity (CondOpt2)
  * Custom rank-one spectral density (CustomOne) -- the energy-resolved density $\mathrm{Tr}[A\,\delta(E-H)]$ of a user-defined `#!python kite.custom.Vertex` operator $A$ from [`#!python calculation.custom_one()`][calculation-custom_one]
  * Custom rank-two Kubo-Bastin response (CustomTwo) -- the Fermi-energy-resolved response reconstructed from the double-Chebyshev moment matrix $\Gamma_{mn}=\mathrm{Tr}[T_m(\hat H)\,B\,T_n(\hat H)\,A]$ of a user-defined vertex pair $A,B$ from [`#!python calculation.custom_two()`][calculation-custom_two]

The SingleShot DC conductivity does not require the post-processing through KITE-tools.


### Advanced usage

KITE-tools supports a set of command-line instructions to force it to use user-specified parameters for each of the quantities mentioned in the previous section. The syntax is as following:

``` bash
./build/KITE-tools archive.h5 --quantity_to_compute1 -key_1 value_1 -key_2 value_2 --quantity_to_compute2 -key_3 value_3 ...
```

Each function to compute is specified after the double hyphens — and the parameters of each function is specified after the single hyphen -. The list of available commands is as follows:

| Function            | Parameter    | Description                                                                                         |
|---------------------|--------------|-----------------------------------------------------------------------------------------------------|
| `#!bash --LDOS`     | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --LDOS`     | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --LDOS`     | `#!bash -K`  | Kernel to use (jackson/green). green requires broadening parameter. Example: `#!bash -K green 0.01` |
| `#!bash --LDOS`     | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --ARPES`    | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --ARPES`    | `#!bash -E`  | min max num Number of energy points                                                                 |
| `#!bash --ARPES`    | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --ARPES`    | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --ARPES`    | `#!bash -V`  | Wave vector of the incident wave                                                                    |
| `#!bash --ARPES`    | `#!bash -O`  | Frequency of the incident wave                                                                      |
| `#!bash --ARPES`    | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --DOS`      | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --DOS`      | `#!bash -E`  | Number of energy points                                                                             |
| `#!bash --DOS`      | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --DOS`      | `#!bash -K`  | Kernel to use (jackson/green). green requires broadening parameter. Example: `#!bash -K green 0.01` |
| `#!bash --DOS`      | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CondDC`   | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondDC`   | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondDC`   | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondDC`   | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondDC`   | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondDC`   | `#!bash -d`  | Broadening parameter of the Dirac delta                                                             |
| `#!bash --CondDC`   | `#!bash -F`  | min max numRange of Fermi energies. min and max may be omitted if only one is required              |
| `#!bash --CondDC`   | `#!bash -t`  | Number of threads                                                                                   |
| `#!bash --CondDC`   | `#!bash -I`  | If `#!bash 0`, CondDC uses the DOS to estimate the integration range                                |
| `#!bash --CondDC`   | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CondOpt`  | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondOpt`  | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondOpt`  | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondOpt`  | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondOpt`  | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --CondOpt`  | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondOpt`  | `#!bash -O`  | min max num Range of frequencies                                                                    |
| `#!bash --CondOpt`  | `#!bash -C`  | `num_d num_g` Recompute the conductivity using only `num_d` Dirac-delta and `num_g` Green's-function moment blocks, for convergence testing (see [Ground Rules][ground_rules]) |
| `#!bash --CondOpt`  | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CondOpt2` | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondOpt2` | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondOpt2` | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondOpt2` | `#!bash -R`  | Ratio $r$ of the second frequency to the first one, $\omega_2=r\,\omega_1$. Default `#!bash 1.0` (second-harmonic/sum-frequency generation). **Setting `#!bash -R -1` switches to the degenerate $\omega_2=-\omega_1$ case, i.e. the DC photocurrent/shift-current (photoconductivity) response**, using a different internal calculation path (`calculate_photo()` instead of `calculate_general()`) |
| `#!bash --CondOpt2` | `#!bash -P`  | If set to 1: writes all the different contributions to separate files                               |
| `#!bash --CondOpt2` | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondOpt2` | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --CondOpt2` | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondOpt2` | `#!bash -O`  | min max num Range of frequencies                                                                    |
| `#!bash --CondOpt2` | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CustomOne` | `#!bash -N` | Name of the output file                                                                            |
| `#!bash --CustomOne` | `#!bash -E` | min max num Number of energy points                                                                |
| `#!bash --CustomOne` | `#!bash -X` | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CustomTwo` | `#!bash -N` | Name of the output file                                                                            |
| `#!bash --CustomTwo` | `#!bash -E` | Number of energy points used in the integration                                                     |
| `#!bash --CustomTwo` | `#!bash -T` | Temperature                                                                                          |
| `#!bash --CustomTwo` | `#!bash -S` | Broadening parameter of the Green's function                                                        |
| `#!bash --CustomTwo` | `#!bash -d` | Broadening parameter of the Dirac delta                                                             |
| `#!bash --CustomTwo` | `#!bash -F` | min max num Range of Fermi energies. min and max may be omitted if only one is required              |
| `#!bash --CustomTwo` | `#!bash -X` | Exclusive. Only calculate this quantity                                                             |

All the values specified in this way are assumed to be in the same units as the ones used in the configuration file. All quantities are double-precision numbers except for the ones representing integers, such as the number of points. This list may be found in KITE-tools, run `KITE-tools --help`:

## Output

In the table below, we specify the name of the files that are created by KITE-tools according to the calculated quantity and the format of the data file.

| Quantity                          | File                         | Column 1          | Column 2         | Column 3         |
|-----------------------------------|------------------------------|-------------------|------------------|------------------|
| Local Density of States           | `#!bash ldos{E}.dat`         | lattice position  | LDOS ($Re$)      |                  |
| ARPES                             | `#!bash arpes.dat`           | k-vector          | ARPES ($Re$)     |                  |	
| Density of States                 | `#!bash dos.dat`             | energy            | DOS ($Re$)       | DOS ($Im$)       |
| Optical Conductivity              | `#!bash optical_cond.dat`    | Frequency         | Opt. Cond ($Re$) | Opt. Cond ($Im$) |
| DC Conductivity                   | `#!bash condDC.dat`          | Fermi energy      | Cond ($Re$)      | Cond ($Im$)      |
| Second-order optical conductivity | `#!bash nonlinear_cond.dat`  | Frequency         | NL Cond ($Re$)   | NL Cond ($Im$)   |
| Single-shot DC Conductivity       | `#!bash [HDF5-filename].dat` | Fermi energy      | Cond ($Re$)      | Cond ($Im$)      |
| Custom rank-one spectral density   | `#!bash custom_one.dat`      | energy            | $\mathrm{Tr}[A\,\delta(E-H)]$ ($Re$) |      |
| Custom rank-two Kubo-Bastin response | `#!bash custom_two.dat`    | Fermi energy      | response ($Re$)  |                  |

* All linear conductivities are in units of $e^2/h$
* Both Planck’s constant and electron charge are set to 1.
* LDOS outputs one file for each requested energy. The energy is in the E in the file name.

For more details on the type of calculations performed during post-processing, check [Resources][resources] where we discuss our method.

<span id="kite-tools-customone"></span>

!!! info "`#!bash --CustomOne`: why the result can be real even if the vertex "looks" imaginary"

    `#!python calculation.custom_one()` builds its vertex $A$ out of building blocks including KITE's
    raw `#!python "vx"`/`#!python "vy"` velocity token, which is $[\hat H,\hat r]$ -- **missing** the
    $i/\hbar$ factor of the textbook Hermitian velocity operator (this is a deliberate KITE convention:
    it keeps the underlying arithmetic simpler, and the compensating factor of $i$ is resolved once,
    here in post-processing, rather than by every caller separately). Whether $\mathrm{Tr}[T_n(\hat H)A]$
    comes out real or purely imaginary depends on the **parity** of how many `#!python "v"`-type tokens
    $A$ is built from — even, and it's real; odd, and it's genuinely imaginary, not an error. KITE
    auto-detects this count directly from the vertex definition (`#!python calculation.custom_one()`
    counts `#!python "v"` tokens itself and stores it as `#!python NumVelocities` in the output file), and
    `#!bash --CustomOne` uses it to rotate an odd-parity result back onto the real axis before applying
    the usual Jackson-kernel reconstruction -- you do not need to insert any compensating `#!python i`
    into the vertex yourself. See [Custom Vertex Operators][custom-vertex-example] for the general
    mechanism, and `#!python examples/haldane_orbital_magnetization.py` for a worked example with an odd
    (single) velocity operator.

<span id="kite-tools-customtwo"></span>

!!! info "`#!bash --CustomTwo`: a general reconstruction, not a spin-Hall-specific formula"

    `#!python calculation.custom_two()` writes a raw double-Chebyshev moment matrix
    $\Gamma_{mn}=\mathrm{Tr}[T_m(\hat H)\,B\,T_n(\hat H)\,A]$ for an arbitrary vertex pair
    $A,B$. As with `#!bash --CustomOne`, each raw stored operator is $X=i^{n_X}\hat X_H$ with
    $\hat X_H$ Hermitian and $n_X$ counting KITE's "missing factor of $i$" per velocity token, so
    $\Gamma_{mn}$ picks up an overall $i^p$ with $p=$ `#!python NumVelocities` (the combined count
    from both $A$ and $B$). Since $i^p$ has **period 4, not 2**, which component of the Chebyshev
    contraction $Z(E)=\sum_{mn}\delta_m(E)\,\Gamma_{mn}\,\mathrm{dgreen}_n(E)$ is physical depends
    on $p\bmod4$:

    | $p\bmod4$ | correct component |
    |---|---|
    | 0 | $-2\,\mathrm{Im}[Z]$ |
    | 1 | $+2\,\mathrm{Re}[Z]$ |
    | 2 | $+2\,\mathrm{Im}[Z]$ |
    | 3 | $-2\,\mathrm{Re}[Z]$ |

    `#!bash --CustomTwo` reads `#!python NumVelocities` directly from the file and applies this
    table automatically -- no manual branch selection needed. This is the same rule used by
    [`#!python examples/rashba_edelstein_graphene_process.py`][ree-example]'s Python
    post-processing (for a $p=1$ vertex pair) and by
    [`#!python examples/kane_mele_spin_hall_process.py`][custom-vertex-example] (for the $p=2$
    spin-Hall pair, where it reduces to the previously-published formula). See
    [Rashba-Edelstein Effect][ree-example] for the full derivation and a worked example.

    **Units**: for `#!python NumVelocities=2` (e.g. spin Hall), the result is in units of
    $e^2/h$, same as `#!bash --CondDC`. For `#!python NumVelocities != 2`, an additional factor
    of `#!bash EnergyScale`$^{N_v-2}$ is applied automatically.

!!! info "Processing the single-shot DC conductivity"

    The single shot DC conductivity does not need any post-processing as it is an energy dependent calculation where the conductivity is calculated on the fly.
    In this particular case, the data is extracted directly from the hdf file with the following python script in the `#!bash kite/tools/`-folder:
    
    ``` bash
    python3 process_single_shot.py output.h5
    ```
    
    The output of this script will be a data-file with name `#!bash output.dat` or similar, structured as given in the table above.

## Examples


### Example 1

``` bash
./KITE-tools h5_file.h5 --DOS -E 1024
```

Processes the .h5 file as usual but ignores the number of energy points in the density of states present there. Instead, KITE-tools will use the value 1024 as specified in the example.

### Example 2

``` bash
./KITE-tools h5_file.h5 --CondDC -E 552 -S 0.01
```

Processes the .h5 file but uses 552 points for the energy integration and a broadening parameter of 0.01.

### Example 3

``` bash
./KITE-tools h5_file.h5 --CondDC -T 0.4 -F 500
```

Calculates the DC conductivity using a temperature of 0.4 and 500 equidistant Fermi energies spanning the spectrum of the Hamiltonian.

### Example 4

``` bash
./KITE-tools h5_file.h --CondDC -F -1.2 2.5 30 --CondOpt -T 93
```

Calculates the DC conductivity using 30 equidistant Fermi energies in the range `#!python [-1.2, 2.5]` and the optical conductivity using a temperature of 93.

### Example 5

``` bash
./KITE-tools h5_file.h5 --CustomOne -E -3 3 400
```

Reconstructs the `#!python custom_one()` spectral density $\mathrm{Tr}[A\,\delta(E-H)]$ over 400 energy points in the range `#!python [-3, 3]`.

### Example 6

``` bash
./KITE-tools h5_file.h5 --CondOpt2 -R -1 -P 1 -O 0 2 200
```

Computes the **DC photocurrent / shift-current response** (the degenerate $\omega_2=-\omega_1$ limit of the
second-order optical conductivity, selected via `#!bash -R -1`) over the frequency range `#!python [0, 2]`
with 200 points, writing each Gamma-contribution to its own file (`#!bash -P 1`).

### Example 7

``` bash
./KITE-tools h5_file.h5 --CustomTwo -F -2 2 200 -T 0.01
```

Reconstructs the `#!python custom_two()` response as a function of Fermi energy over 200
equidistant points in `#!python [-2, 2]`, at temperature `#!bash 0.01`.

@@include[kite_tools_readme.md](kite_tools_readme.md)

[resources]: ../background/index.md
[ground_rules]: ../documentation/optimization.md
[calculation-custom_one]: kite.md#calculation-custom_one
[calculation-custom_two]: kite.md#calculation-custom_two
[custom-vertex-example]: ../documentation/examples/custom_vertex_operators.md
[ree-example]: ../documentation/examples/rashba_edelstein.md
 