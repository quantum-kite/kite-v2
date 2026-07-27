# What changed and why (2026-07-26 hopping/phase-convention fix)

Pushed to `quantum-kite/master` as two commits: `7440dd5` and `3628387`.
This file is a plain-language index of every file touched, why, and — just
as important — what was **not** touched, since `build_planewave()` (the
ARPES plane-wave code) is discussed at length below but was never edited.

## Where this started

An external finding from a sibling project reported that
`Lattice.add_one_hopping`/`add_hoppings` store `H[row=from_sub, col=to_sub]`
— i.e. `from_sub` indexes the row and `to_sub` the column, the opposite of
how "hop from A to B" reads intuitively — and that this was completely
undocumented. You asked me to review every complex-hopping example against
this convention and check the docs actually match the code. That review
surfaced a real bug, which a follow-up third-party audit then found I had
only half-fixed. Both rounds are covered below.

## Files touched, and why

### `src/kite/lattice.py`
**Docstrings only, no logic change.** Added an explicit statement of the
row/col convention to `add_one_hopping`/`add_hoppings` (this was the
original ask — the convention was real and correct, just unwritten).

### `src/kite/visualize.py`
**The actual bug fix.** `hamiltonian_k()` is a pure-Python preview/plotting
helper used to draw a band structure from a `Lattice` object directly in
Python — it is a separate code path from the real simulation, which goes
through `config_system()` → HDF5 → the C++ engine. `hamiltonian_k()` had
two independent bugs:
1. It built `H0[to_id, from_id]` instead of `H0[from_id, to_id]` — backwards
   relative to `config_system()`'s real export.
2. It used `exp(-i k·r)` instead of `exp(+i k·r)` — the wrong sign of k,
   i.e. it silently computed `H(-k)` instead of `H(k)`.

I found and fixed (1) myself while doing the review you asked for. I fixed
(2) after the external audit caught it — I had missed it because both bugs
are invisible to any check that only looks at eigenvalues (band energies,
gaps): `H` and `H^T` share eigenvalues, and `H(k)`/`H(-k)` share eigenvalues
too, so every "does this reproduce the known band structure" check I ran
passed regardless. To confirm the correct sign, I read KITE's actual C++
ARPES plane-wave state, `build_planewave()` in
`Src/Vector/KPM_Vector2D.cpp`/`KPM_Vector3D.cpp` — **I only read this file
as a ground-truth reference, I did not modify it.** It uses `exp(+i k·r)`,
which is what `hamiltonian_k()` now matches.

**Nothing about real ARPES simulations changed.** `hamiltonian_k()` is not
used by the ARPES example or by `KITEx`/`KITE-tools`; this bug never
affected any real simulation output, only this one Python plotting helper.

### `tests/hamiltonian_k_convention_regression.py`
**Rewritten test.** My first version of this test only used a "paired"
hopping fixture (both A→B and B→A stored explicitly, like the Weyl model),
which structurally cannot distinguish `exp(+ik)` from `exp(-ik)` — the
phase collapses to a cosine either way. That's why it passed even with bug
(2) still present. Added a second, "unpaired" fixture (one direction
stored, the mirror auto-generated) that actually exercises the phase sign,
matching the audit's minimal reproduction.

### `examples/process_haldane_orbital_magnetization.py`
**Physics bug fix, not a documentation change.** This script independently
reimplements a Bloch Hamiltonian (for a Berry-curvature cross-check) using
the same row/col and phase convention as `hamiltonian_k()`, so it carried
its own copies of both bugs. It also had a separate, unrelated sign error in
its modern-theory orbital-magnetization formula (`+Im[...]` instead of
`-Im[...]`). The two errors happened to cancel, so the script's output
matched the documented result (a Chern number of +1) for the wrong reason.
Fixing both together changes the correct physical answer to **C = -1**
(and the slope relation from `C/(2π)` to `-C/(2π)`), while reproducing the
same numerical curve — I reran the real `KITEx`/`KITE-tools` pipeline
end-to-end afterward to confirm this (slope 0.1574 vs. predicted 0.1592,
1.1% agreement).

### `examples/haldane_orbital_magnetization.py`
**One-line docstring fix.** The vertex-operator derivation had dropped an
`i` when simplifying (`... = -(e/(2ħc·Area))·A` should be `·i·A`, since `A`
is anti-Hermitian and the `i` isn't optional).

### `examples/README.md` and `docs/documentation/examples/orbital_magnetization.md`
**Doc text updated to match the corrected C=-1 result** described above —
no code in either file.

### `docs/api/kite.md`
**Two spots fixed:** the "gauge convention" description of `H(k)` (updated
to the correct row/col + `exp(+ik)` form, with the derivation spelled out),
and the Weyl-semimetal worked example's surrounding text, which had
incorrectly called a sign flip "just a basis choice" — it's actually a
physical change in the Weyl node's chirality. **This page's ARPES section
was not touched.**

### `maintenance/native-lattice-viz-plan.md`
**Root-cause fix.** This is the original design note for the plotting
feature, written before any code existed. It already had both the row/col
label and the phase sign backwards from the start — that's where both bugs
actually originated. Corrected to match the now-confirmed convention.

### `maintenance/2026-07-26-hopping-convention-audit.md`
**My own running log of this investigation**, updated to record the
follow-up audit's findings as resolved (this is an internal maintenance
note, not user-facing documentation).

## What was *not* touched

- **ARPES**: `examples/arpes_tmd.py`, `docs/documentation/examples/spectral_function.md`,
  and `build_planewave()` itself (`Src/Vector/KPM_Vector2D.cpp`/`KPM_Vector3D.cpp`)
  were all read-only references used to confirm the correct sign convention —
  none were edited. ARPES simulation output is unaffected by any of this.
- Any other complex-hopping example (Rashba-Edelstein, Fu-Kane-Mele, Kane-Mele,
  AFM ribbon, the other Weyl scripts) — checked against the confirmed convention
  and found already correct; no changes needed.
- `config_system()` and the real HDF5-export path — always correct; the bug was
  confined to the separate `hamiltonian_k()` Python preview helper and the one
  script that copied its (buggy) convention.
