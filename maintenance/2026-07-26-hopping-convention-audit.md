# Hopping row/col convention audit (2026-07-26)

Triggered by an external finding shared from a sibling project (t2g-orbital-relaxation's
`KITE_HOPPING_CONVENTION_FINDING.md`): `Lattice.add_one_hopping`/`add_hoppings` set
`value = H[row=from_sub, col=to_sub]` — `from_sub` indexes the ROW, `to_sub` the COLUMN,
the opposite of the "hop from A to B" reading some conventions use — and this was
completely undocumented. Confirmed directly against `config_system()`'s actual sparse
matrix construction before touching anything.

## Fixed

- **`src/kite/lattice.py`**: `add_one_hopping`/`add_hoppings` docstrings now state the
  row/col convention explicitly (the external finding's suggested fix).
- **A real, independently-discovered bug in `kite.visualize.hamiltonian_k()`**: it built
  `H0[to_id, from_id]` — backwards relative to `config_system()`. Confirmed three ways:
  hand-reconstruction of `config_system`'s real-space export (exact transpose match),
  tracing the actual C++ propagation code (`Src/Vector/KPM_Vector2D.cpp:520`,
  `phi0[io] += hopping(ib,io)*phiM1[dist(ib,io)]`, i.e. row=`from_id`), and an exact
  elementwise match after fixing it. **Never affected real KITEx simulation results** —
  `config_system()` builds the real HDF5 export independently and was always correct.
  Only affected this preview/verification helper and anything that duplicated its
  convention instead of comparing against real KITEx output.
- Added `tests/hamiltonian_k_convention_regression.py` and wired it into CI
  (`build-and-test`, pure Python, no KITEx needed) — this bug is invisible to any
  eigenvalue-based check (H and H^T share eigenvalues), which is exactly how it went
  undetected through every prior "matches known closed form" verification. The test
  compares matrix elements against an independent reconstruction, not eigenvalues.
- Fixed `docs/api/kite.md`'s two copies of the stale convention description and the
  Weyl-semimetal closed-form claim, which changes from
  `H(k)=t[cos(kx)sx+cos(ky)sy+cos(kz)sz]` to `t[cos(kx)sx-cos(ky)sy+cos(kz)sz]` — the
  same physical model up to relabeling `sy -> -sy`, not a different system, but the
  literal formula stated no longer matched what `hamiltonian_k` (now fixed) computes.

## Reviewed, no bug found

- **`rashba_edelstein_graphene.py`**: derived the correct Rashba SOC matrix element by
  hand under the confirmed convention; matches the code's `up_to_dn`/`dn_to_up`
  functions exactly (up to a real scale factor).
- **`weyl_lt.py`, `dos_t_symmetric_cubic_weyl_sm.py`, `optcond_t_symmetric_cubic_weyl_sm.py`**:
  exact match to the (corrected) closed form after the `hamiltonian_k` fix.
- **`dos_fu_kane_mele_model.py`**: exact Kramers-pair degeneracy confirmed numerically
  at multiple k-points (a structural check independent of eigenvalue-only checks).
- **Haldane, Kane-Mele, AFM ribbon same-sublattice complex hoppings**: structurally
  immune — `from_sub == to_sub` makes the row/col swap a no-op on a diagonal entry.
- `kite.repository`'s TMD builder: no complex hoppings at all in the current model.

## Open — not resolved, flagged for follow-up

`examples/process_haldane_orbital_magnetization.py` has its own independent copy of the
same (now-fixed-elsewhere) convention, used to cross-check the real, actually-computed
KITEx `custom_one` orbital-magnetization output against a topological Chern-number
prediction (documented as C=+1). Ran the actual KITEx computation: the real data gives a
magnetization slope of ~+0.156 in the gap, matching C=+1 as documented. Recomputing the
Chern number via a Fukui-Hatsugai-Suzuki lattice calculation on the now-fixed
`hamiltonian_k` gives C=-1 instead — a genuine contradiction. Likely explanation: the
`custom_one` vertex `A = x*H*y - y*H*x` (built in the C++ code) has its own
row/col-sensitive behavior that doesn't simply compose with `hamiltonian_k`'s convention
the way assumed. Did not modify this file's convention without resolving the
contradiction first — spawned as a follow-up task (`task_e72ed758`) rather than guessing.
