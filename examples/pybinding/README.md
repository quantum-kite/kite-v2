# Pybinding-based examples

Every script in this folder builds its lattice with [Pybinding](https://docs.pybinding.site/),
which is an optional extra for KITE (`pip install -e ".[pybinding]"`), not a core dependency.

For most of these, the *only* thing pybinding is used for is lattice construction via `pb.Lattice(...)`,
`.add_sublattices(...)`, and `.add_hoppings(...)` — KITE's own `kite.lattice.Lattice` class
(`src/kite/lattice.py`) mirrors that same API and needs no extra dependency. `weyl_pb.py` in this folder
and `../weyl_lt.py` at the top level are the same physical system, one written each way, if you want to
compare them directly.

A few examples here use more of pybinding than just `pb.Lattice`:

- `dos_mixed_disorder.py`, `dos_on_site_disorder.py`, `dos_optcond_gaussian_disorder.py` use
  `pybinding.repository.graphene`, a built-in lattice template, instead of defining the lattice by hand.
- `dos_twisted_bilayer.py` uses pybinding's own `Model`/`kpm`/`translational_symmetry` system directly —
  a fundamentally different workflow from KITE's, not just lattice construction.

`arpes_bilayer.py`, `arpes_cubic.py`, and `arpes_tmd.py` used to be in this folder too, for
`pb.results.make_path(...)`, a k-path-interpolation helper — now reproduced natively as
`kite.visualize.make_path(...)` (`src/kite/visualize.py`), so all three have moved to the top level
of `examples/`. `arpes_tmd.py`'s locally-defined `tmd_monolayer()` was also replaced by
`kite.repository.group6_tmd.monolayer_3band()` (`src/kite/repository.py`).

A native `kite.repository.graphene.monolayer(...)` template also exists now (verified against
`pybinding.repository.graphene`) but hasn't yet been wired into the three disorder examples above —
that conversion, and deciding what (if anything) to do about `dos_twisted_bilayer.py`'s fundamentally
different workflow, are what's left to shrink this folder further.
See `../../maintenance/installability-plan.md` for the broader context on why pybinding was made optional.
