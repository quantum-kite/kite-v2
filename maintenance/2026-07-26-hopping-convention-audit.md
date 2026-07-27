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

## Resolved (2026-07-26, follow-up pass): the phase sign and the Haldane contradiction

The "open" item below (as originally written in this file) was itself based on an
incomplete fix: `hamiltonian_k()`'s row/col was corrected, but a SEPARATE bug —
the phase sign, `exp(-i k.r)` instead of `exp(+i k.r)` — was still present, and my own
speculation that the `custom_one` vertex had a hidden row/col-sensitive behavior was
wrong. An independent audit found and resolved both:

- **Phase sign**: `hamiltonian_k()`'s `exp(-1j*...)` computed `H(-k)` instead of `H(k)`.
  Confirmed directly against KITE's actual C++ ARPES plane-wave state
  (`build_planewave()` in `Src/Vector/KPM_Vector2D.cpp`/`KPM_Vector3D.cpp`: `exp(+i
  k.(r+d))`, a plus sign, not assumed). Fixed to `exp(+1j*...)`. Like the row/col bug,
  this is invisible to any eigenvalue-only check (H(k) and H(-k) share a spectrum at any
  k where -k is also sampled) and to the row-fix's own regression test, whose Weyl
  fixture used PAIRED hoppings (both directions stored) whose phase dependence collapses
  to cosines and can't distinguish +ik from -ik. Added an UNPAIRED complex-hopping test
  case to `tests/hamiltonian_k_convention_regression.py` to actually catch this.
- **No hidden transpose in `custom_one`**: confirmed directly against
  `Src/Hamiltonian/HamiltonianRegular.cpp:127-136` — the velocity vertex `v1(ih,io)` uses
  the exact same row=`from_id` indexing as the hopping matrix itself. My earlier
  speculation about this being the culprit was wrong.
- **`process_haldane_orbital_magnetization.py`** had its own independent copy of BOTH the
  original row/col bug and (implicitly, via the same `hamiltonian_k`-derived formula) the
  phase-sign issue, used for its Berry-curvature cross-check — fixed to match the
  corrected convention (row=`from_id`, col=`to_id`, `exp(+i k.r)`). Separately, its
  modern-theory formula had the opposite-of-correct overall sign for the physical
  electron-charge magnetization (`+Im[...]` where the correct expression is `-Im[...]`).
  These two errors were compensating: fixing only one would have visibly broken the
  documented agreement with the real KITEx data; fixing both together reproduces the
  existing curve exactly (RMS difference to the pre-fix curve: negligible, well below
  plotting resolution) while correcting the physical interpretation from **C=+1 to
  C=-1** and the slope relation from **C/(2π) to -C/(2π)**. Verified end-to-end: ran the
  real KITEx computation fresh, confirmed the actual data (slope ≈ +0.157) matches
  `-C/(2π) = +0.159` for `C=-1` to ~1%. Updated `CHERN_NUMBER`, `EXPECTED_SLOPE`, the plot
  title, and every doc/README reference to this example's Chern number and slope
  relation accordingly.
- Also fixed: `haldane_orbital_magnetization.py`'s own docstring dropped an `i` when
  rewriting the operator prefactor (`M_z = -(e/(2hbar*c*Area))*A` instead of `*i*A` — A is
  anti-Hermitian, so the `i` cannot be dropped); and `maintenance/native-lattice-viz-plan.md`,
  the original source of both the row/col and phase-sign conventions this bug traces back
  to, which had them backwards from the start.

This closes the item without needing to guess: every claim above was verified directly
against the actual C++ source (not re-derived from scratch) before any code changed.
