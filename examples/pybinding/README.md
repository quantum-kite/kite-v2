# Pybinding-based examples

Every script in this folder builds its lattice with [Pybinding](https://docs.pybinding.site/),
which is an optional extra for KITE (`pip install -e ".[pybinding]"`), not a core dependency.

For most of these, the *only* thing pybinding is used for is lattice construction via `pb.Lattice(...)`,
`.add_sublattices(...)`, and `.add_hoppings(...)` — KITE's own `kite.lattice.Lattice` class
(`src/kite/lattice.py`) mirrors that same API and needs no extra dependency. `weyl_pb.py` in this folder
and `../weyl_lt.py` at the top level are the same physical system, one written each way, if you want to
compare them directly.

A few examples here use more of pybinding than just `pb.Lattice`:

- `arpes_bilayer.py`, `arpes_cubic.py`, `arpes_tmd.py` also use `pb.results.make_path(...)`, a small
  k-path-interpolation helper.
- `dos_mixed_disorder.py`, `dos_on_site_disorder.py`, `dos_optcond_gaussian_disorder.py` use
  `pybinding.repository.graphene`, a built-in lattice template, instead of defining the lattice by hand.
- `dos_twisted_bilayer.py` uses pybinding's own `Model`/`kpm`/`translational_symmetry` system directly —
  a fundamentally different workflow from KITE's, not just lattice construction.

We plan to reproduce the two or three specific pybinding functions these examples actually need
(the k-path helper and the graphene lattice template being the main ones) directly in
`kite.lattice`/`kite.custom`, so the rest of these examples can eventually move out of this folder too.
See `../../maintenance/installability-plan.md` for the broader context on why pybinding was made optional.
