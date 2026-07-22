!!! Info

    This file describes the contents of the
    [`#!bash kite/examples/custom_local_potential/`][examples-clp-github]-folder.

This example shows how to use KITE's custom local potential functionality: a way to add an arbitrary,
deterministic, site- and orbital-dependent on-site energy to the tight-binding Hamiltonian — for example a
gate-defined electrostatic confinement potential, a p-n junction profile, or a staggered/orbital-selective
mass term — on top of (or instead of) the stochastic disorder provided by [`#!python kite.Disorder`][disorder].

!!! warning "This mechanism is compiled into KITEx, not loaded at runtime"
    An older version of this example described compiling a standalone `custom_potential.cpp` into a shared
    library `libaux.so.1` that KITEx would `dlopen` at runtime. **That mechanism no longer exists in the
    current codebase**: there is no `dlopen`/`dlsym`/shared-library loading anywhere in KITEx, and
    `custom_potential.cpp` in this folder `#!cpp #include "aux.hpp"`, a header that does not exist anywhere
    in the current `Src/` tree — so following the old instructions would fail to compile. The
    `custom_potential.cpp`, `custom_potential_default.cpp`, `test.sh`, `test_potential.py` and
    `test_potential.ipynb` files still present in this folder are leftovers from that old approach; they use
    the old, now-nonexistent `#!cpp double uncorrelated_wrapper(double *vec, unsigned orb)` free-function
    API and `test.sh` copies its file into `Src/Hamiltonian/custom_potential.cpp`, a path that has never
    matched where the real implementation lives (see below). Do not use them as a template — use
    `Src/custom.cpp`, described below, instead.

## How it actually works today

The potential is implemented as a single C++ member function that you edit directly and recompile:

``` cpp
// Src/custom.cpp
template <typename T, unsigned D>
void Hamiltonian<T, D>::assign_local_potential(
  const unsigned index_,
  const unsigned io_
)
```

* `index_` is a **local** (ghost-padded, sub-domain) lattice-site index. Inside the function it is converted
  to a global unit-cell coordinate (`r.convertCoordinates(global, local)`) and then to a Cartesian position
  via `#!cpp pos = r.rLat * v_global` (`r.rLat` holds the lattice vectors as columns).
* `io_` is the orbital index, so the function can also branch on which orbital is being evaluated.
* The function is expected to set the member `custom_pot` (of type `value_type`, i.e. the same
  floating-point/real type as the rest of the Hamiltonian) as a side effect. `custom_pot` is expressed in
  the **same physical energy units as the tight-binding hoppings** — the conversion to KITE's internal
  rescaled units happens automatically afterwards.

!!! warning "`pos` is the unit-cell origin, not the orbital's own position"
    `pos = r.rLat * v_global` only gives the position of the **origin of the unit cell**, not the actual
    position of orbital `io_` within that cell. For a lattice where different orbitals of the same
    sublattice/unit-cell sit at different physical locations, a potential that is meant to depend on the true
    intra-cell position needs to add the orbital's own offset, `#!cpp r.rOrb.col(io_)`, to `pos` itself. If
    you only need orbital-dependent behavior (not sub-unit-cell spatial resolution), branching on `io_`
    alone, as in the example below, is enough and this correction isn't needed.

`Src/custom.cpp` is compiled directly into the `KITEx` executable by CMake's recursive glob
(`FILE(GLOB_RECURSE SRCFILES Src/*.cpp)`, `CMakeLists.txt`), together with every other `.cpp` file under
`Src/`. There is no separate plugin step: to use your own potential you edit `Src/custom.cpp` in place (or
back up the original and drop in a replacement with the same function signature) and do a normal full
rebuild:

``` bash
cd build
cmake .. && make
```

## The shipped example: a soft-wall (gate-defined) quantum dot

The `assign_local_potential` currently checked in to `Src/custom.cpp` implements a smooth, radially
symmetric confinement potential — physically, the kind of "soft quantum dot" profile produced by an
electrostatic gate: zero potential inside a disk around the center of the sample, rising smoothly (rather
than as a hard wall) outside it. For a 2D system (`D == 2`):

``` cpp
const Eigen::Matrix<double, D, 1> origin(0.5 * r.Lt[0], 0.5 * r.Lt[1]);
const Eigen::Matrix<double, D, 1> rel_pos = pos - r.rLat * origin;
const value_type radius = rel_pos.norm();
const value_type boundary = r.Lt[0] / 8;
if (radius < boundary)
  custom_pot = 0.0;
else {
  const value_type sigma = r.Lt[0] / 8;
  const value_type Vmax = 2.0;
  const value_type dr = radius - boundary;
  custom_pot = Vmax * (1.0 - std::exp(-0.5 * dr * dr / (sigma * sigma)));
}
```

`r.Lt[0]`/`r.Lt[1]` are the total sample dimensions (in unit cells) along each direction, so `origin` is the
sample's geometric center. `boundary = r.Lt[0] / 8` sets the radius of the field-free disk to one eighth of
the sample's linear size; inside that disk `custom_pot` is exactly zero. Outside it, `custom_pot` rises
smoothly from zero towards an asymptotic value `Vmax = 2.0` (in the same energy units as the hoppings) with
a Gaussian-CDF-like envelope of width `sigma = r.Lt[0] / 8` — so the potential is essentially "off" at the
boundary and reaches within `1/e`-type closeness of `Vmax` a distance `~sigma` further out, avoiding the
sharp features (and the associated numerical/Gibbs artifacts) a hard step would introduce. There is no
`io_`-dependence in the current example — every orbital sees the same confinement — and no `D == 3` branch
is implemented (the function is a no-op for 3D lattices unless you add one).

Other physically meaningful uses of this same mechanism, which you can implement by writing a different
body for `assign_local_potential` (not provided as ready-made examples here):

* **p-n junctions**: a step- or tanh-shaped potential along one Cartesian direction of `pos`, instead of a
  radial profile.
* **Orbital-selective / staggered mass terms**: branch on `io_` to give different sublattices or orbitals
  different constant (or position-dependent) on-site energies.
* **Deterministic disorder landscapes**: a spatially-varying but non-random potential, as an alternative (or
  complement) to the stochastic Gaussian/Uniform disorder from `kite.Disorder`.

## Combining with `kite.Disorder`

The custom potential is wired in at `Src/Hamiltonian/Hamiltonian.cpp` (`distribute_AndersonDisorder()`, in
the `#!cpp if (check_local_potential) { ... }` block): for every site and orbital,
`#!cpp custom_pot / EnergyScale` is **added directly on top of** whatever on-site Anderson term
(`U_Anderson`) is already there from a `kite.Disorder` configuration. In other words this is additive, not a
replacement: you can define a deterministic spatial potential profile with this mechanism *and* superimpose
stochastic Gaussian/Uniform on-site disorder from `kite.Disorder` in the same calculation, and both
contributions will be present simultaneously in the final on-site energies.

## Python-facing configuration

`config.py` in this folder builds the rest of the `.h5` file the usual way, via `kite.Configuration(...)`
and `kite.config_system(...)`. The flag that actually turns the mechanism on in KITEx is
`#!python custom_potential` (an integer, default `0`):

``` python
configuration = kite.Configuration(
    ...,
    custom_potential=1,
)
```

`custom_potential` is written to the HDF5 file as `/Hamiltonian/CustomLocalPotential`, which is what
`Src/Hamiltonian/Hamiltonian.cpp` reads into `check_local_potential` and checks before calling
`assign_local_potential` at all.

!!! warning "`custom_local` / `custom_local_print` are accepted but currently do nothing in KITEx"
    `kite.Configuration` also accepts two older flags, `#!python custom_local` (bool) and
    `#!python custom_local_print` (bool). They are written to HDF5 as `/Hamiltonian/CustomLocalEnergy` and
    `/Hamiltonian/PrintCustomLocalEnergy` respectively — but nothing in the current `Src/` (or `tools/Src/`)
    C++ ever reads either of those two datasets. In particular, `custom_local_print=True` does **not**
    currently make KITEx dump `local_potentialX.dat` files (the old per-thread potential dump this flag used
    to control): the option is accepted for backward compatibility, but produces no output.

    **This previously affected the `config.py` checked in to this folder**: it set `custom_local=True` and
    `custom_local_print=True` but not `custom_potential=1`, so running it would silently *not* activate
    `assign_local_potential` in KITEx (the two "old" flags did nothing, and the one flag that matters
    defaults to off). `config.py` now sets `custom_potential=1` explicitly, so the shipped example does
    activate the soft-quantum-dot potential described above — if you copy this configuration elsewhere,
    remember `custom_potential=1` is the flag that actually matters, not `custom_local`.

## Running the example

``` bash
python config.py
../../build/KITEx config.h5
../../tools/build/KITE-tools config.h5
```

(rebuild KITEx first, `cd build && cmake .. && make`, if you've edited `Src/custom.cpp` to try your own
potential).

[examples-clp-github]: https://github.com/quantum-kite/kite-v2/tree/master/examples/custom_local_potential
[disorder]: ../../api/kite.md#disorder
