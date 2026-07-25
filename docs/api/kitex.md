## Usage

Its default usage is very simple:

``` bash
    ./KITEx archive.h5
```

The program only accepts an HDF5 file as input. No other command-line arguments are accepted.

## Exit behavior

KITEx exits with status `1` (and a one-line message on stdout) in each of the following cases,
checked in this order before any simulation work begins:

- no filename argument was given;
- `/DIM` in the input file is not `1`, `2`, or `3`;
- `/PRECISION` is not `0`, `1`, or `2`;
- `/IS_COMPLEX` is not `0` or `1`;
- the combination of `/PRECISION`, `/DIM`, and `/IS_COMPLEX` selects a real/complex-dimension
  instantiation that this particular build was not compiled with (see
  [Real/complex and precision constraints](#realcomplex-and-precision-constraints) below).

On success it exits with status `0` after printing `Done.`.

## What KITEx reads

KITEx expects the HDF5 file produced by `#!python kite.config_system(...)` (see the
[Python API][kitepython] and [HDF5 editing][editing-hdf5] pages for the full schema). At minimum,
the root of the file must contain:

| Dataset | Meaning |
| --- | --- |
| `/DIM` | Number of spatial dimensions (1, 2, or 3). |
| `/PRECISION` | Floating-point precision: `0` = `float`, `1` = `double`, `2` = `long double`. |
| `/IS_COMPLEX` | `0` for a real Hamiltonian, `1` for complex. |
| `/EnergyScale`, `/EnergyShift` | The affine rescaling `H_tilde = (H - EnergyShift) / EnergyScale` that maps the physical spectrum onto KITE's required Chebyshev domain `[-1, 1]` (with a small safety margin); set from `#!python Configuration`'s `spectrum_range` (or estimated automatically if not given). |
| `/NOrbitals`, `/LattVectors` | Total orbital count and the primitive lattice vectors, used by post-processing to convert results into physical units. |

Beyond that, each requested calculation (`#!python calculation.dos(...)`,
`#!python calculation.ldos_map(...)`, etc.) has its own parameter group under `/Calculation/<name>/`
(for example `/Calculation/dos/NumMoments`), written by `#!python config_system()` and read back by
the corresponding `calc_*()` dispatch routine in `Src/Simulation/GlobalSimulation.cpp`. A group's
absence is how KITEx decides that particular calculation was not requested -- most `calc_*()`
routines gate on a presence check of one of their own parameters (e.g. `NumMoments`) before doing
any work, so leaving a calculation's group out of the file is equivalent to disabling it, not an
error.

## What KITEx writes

Results are written back into the *same* HDF5 file, under `/Calculation/<name>/<ResultDataset>`
(e.g. `/Calculation/dos/lMU` for `dos()`, `/Calculation/ldos_map/Map` for `ldos_map()`) -- KITEx
never creates a separate output file. This is why the three-stage workflow always passes the same
filename to every stage: `#!python kite.config_system(..., filename="x.h5")`, then
`#!bash ./KITEx x.h5`, then `#!bash ./KITE-tools x.h5`.

## Real/complex and precision constraints

KITEx is a template-instantiated C++ program: `Src/main.cpp` selects one of
`float`/`double`/`long double`, each either real or `#!cpp std::complex<...>`, crossed with
1D/2D/3D, purely from `/PRECISION`, `/IS_COMPLEX`, and `/DIM` -- there are 18 possible
instantiations in total. Not all of them are necessarily compiled into a given binary (build flags
can narrow this set to reduce compile time/binary size); requesting a combination that was not
compiled in the target you're running exits with status 1 and the "check if the code has been
compiled with support for complex functions" message. A calculation that has any complex-valued
ingredient (e.g. Peierls/magnetic-field phases, complex hoppings, complex custom operators) must
set `#!python is_complex=True` on `#!python Configuration` -- silently leaving it `False` for such
a model produces incorrect results rather than an error.

## Decomposition

`#!python Configuration`'s `divisions` parameter sets the number of spatial domain-decomposition
parts along each lattice direction; their product is the number of OpenMP threads KITEx will use
(it should not exceed the number of available processor cores). `length` sets the total number of
unit cells along each direction and must be an integer multiple of `divisions` times the internal
tile size the code was compiled with.

## Python / KITEx / KITE-tools compatibility

The three stages must agree on the same HDF5 file and, implicitly, on the same KITE version: the
schema (dataset names, group layout, the rescaling convention above) is a version-coupled contract
between `#!python kite.config_system()` (write), KITEx (read parameters, write raw moments/results),
and KITE-tools (read raw results, reconstruct physical quantities). Mixing a `config_system()` from
one KITE version with a KITEx/KITE-tools binary built from a different version is not supported and
can fail silently (a missing/renamed dataset is generally treated the same as "this calculation
wasn't requested", not as a version-mismatch error) rather than loudly.

[kitepython]: kite.md
[editing-hdf5]: ../documentation/editing_hdf_files.md
