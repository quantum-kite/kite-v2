# Stage-1 third-pass audit corrections checklist

Tracks disposition of every finding in the independent third-pass audit
`kite_v2_stage1_third_pass_audit_2026-07-25.tex`/`.pdf` (audited commit `96bd00f`). Third
pass of external review, following `2026-07-25-stage1-audit-corrections.md` and
`2026-07-25-stage1-reaudit-corrections.md`. Updated as work proceeds; not a substitute
for reading the audit itself.

## Executive findings (audit's own table)

| # | Finding | Severity | Status | Area |
|---|---|---|---|---|
| 1 | Full operator dimension now checked against actual lattice orbital count at export time | P1 | ✅ Already fixed (prior pass) | Core safety |
| 2 | A Hermitian LDOS operator can be mutated into non-Hermitian after the early check and is still exported | P1 | ✅ Fixed | Core safety |
| 3 | Plain and operator LDOS/maps still use separate propagation and disorder/random-vector realizations (documented, not fixed) | P1 | 🟡 Explicit decision recorded: kept as-is, documented more prominently; shared-run redesign is a real follow-up, not attempted | Core semantics |
| 4 | Duplicate orbital indices accepted; silently alias two named orbitals to one matrix entry | P2 | ✅ Fixed | API safety |
| 5 | AFM boundary and broken sublattice-chiral symmetry now described correctly | P2 | ✅ Already fixed (prior pass) | AFM physics |
| 6 | Real-space AFM figure still labels an energy-resolved spectral density as Sz in units of hbar; title says "spin polarization" | P2 | ✅ Fixed | Scientific plotting |
| 7 | Band plot still labels finite-width branch as E=±Delta; prose says only E≈±Delta away from k=0 | P2 | ✅ Fixed | Scientific plotting |
| 8 | CI now runs operator regressions and a strict documentation build | P2 | ✅ Already fixed (prior pass) | CI |
| 9 | New stochastic test has correct numerical expectation but wrong derivation text: writes P_up=(I+sigma_y)/2, should be (I+sigma_z)/2 | P3 | ✅ Fixed | Test documentation |
| 10 | Dimension-test header claims duplicate/non-contiguous cases are exercised; test list contains neither | P3 | ✅ Fixed (coverage extended to match, not just reworded) | Test documentation |
| 11 | Six new commits are narrow, ordered, co-authored, synchronized — strong improvement | P3 | ✅ Already fixed (prior pass) | Organization |

## Gate 1 — finalize operator safety

- [x] Re-run LDOS/LDOS-map label and Hermiticity validation inside `config_system()` on
      the final matrices. Reproduced the audit's adversarial sequence first (register
      Hermitian `l0`, call `ldos_map(operators=['l0'])`, then mutate `l0` non-Hermitian
      via another `add_orbital_coupling()`, then export) — confirmed the bug was real
      (export succeeded with a non-Hermitian operator). Fixed in `config_system()`
      (`src/kite/__init__.py`): collects every label referenced by any `_ldos`/`_ldos_map`
      request and re-checks Hermiticity on the final matrix at export time. Deliberately
      scoped to LDOS-referenced labels only, not every custom operator — `custom_one()`/
      `custom_two()`/`gaussian_wave_packet()` operators can legitimately be non-Hermitian
      intermediates (e.g. velocity operators), so a blanket check would be wrong.
      Re-verified the reproduction is now rejected.
- [x] Validate orbital-index values as unique, contiguous, non-negative integers, cross-
      checked against the lattice orbital count. Reproduced the audit's duplicate-index
      case first (`add_orbital_index('A', 0)` + `add_orbital_index('B', 0)` → couplings
      to A and B silently overwrite the same matrix entry, confirmed the resulting
      matrix was wrong) before fixing. Fixed in `config_system()`: requires the index
      map to be a bijection onto `0..N_orb-1`, raising `ValueError` naming the aliased
      indices/names on duplicates, or the actual vs. expected index set on gaps/negatives.
- [x] Add adversarial regression tests for: duplicate indices, post-registration
      Hermiticity mutation, negative indices, out-of-range indices. Added all four to
      `tests/operator_dimension_regression.py` (now 8 cases total, all passing). Learned
      along the way that non-contiguous/out-of-range indices are usually caught earlier
      by plain numpy bounds-checking inside `add_orbital_coupling()` itself (`IndexError`,
      since the matrix is sized by name-count and an out-of-range index typically exceeds
      that immediately) rather than by the new export-time check — both are legitimate
      rejections, so the test helper accepts either exception type.
- [x] Explicit decision (recorded, not silently deferred): the public API does **not**
      promise plain and operator LDOS share one stochastic realization. Fusing them into
      one shared propagation loop would touch the core KPM dispatch in
      `Src/Simulation/GlobalSimulation.cpp` and both `SimulationLDoS*.cpp` files — a
      much larger, higher-risk change than a documentation fix, and out of proportion to
      what this audit pass should attempt. Kept as two independent calculations, and
      made the cost/randomness distinction more prominent than before: added an explicit
      "known architectural limitation, not a bug" note to both `ldos()`'s and
      `ldos_map()`'s entries in `docs/api/kite.md` (previously only the Python
      docstrings and the AFM example page said this), stating plainly that plain and
      operator maps don't share a realization even within the same call. If this
      limitation is judged unacceptable later, the fix is a genuine redesign, not
      another wording pass — tracked as a real follow-up, not implemented here.

## Gate 2 — finish AFM scientific labeling

- [x] Change the band plot's in-panel annotation from `E=±Delta` to `E≈±Delta`
      (`afm_zigzag_bands.py`'s `ax_ribbon.annotate(...)`); figure regenerated.
- [x] Rename the real-space plotted field/label to `rho_Sz(R,E)` and correct its units
      to `hbar/t` (was `hbar`, wrong for a spectral density). Fixed the colorbar label,
      the vacancy-candidate caption text (`$S_z^A$` -> `$\rho_{S_z}^A$`), and the doc
      page's figure caption in `afm_zigzag_ribbon_process.py` /
      `afm_zigzag_ribbon.md`; figure regenerated and copied into docs assets.
- [x] Replace "edge spin polarization" in the figure title/caption with "spin-weighted
      spectral density" (figure title now: "Spin-weighted spectral density in a
      Néel-gapped bearded/Klein ribbon").
- [x] Compress the off-resonance rationale. Removed the repeated near-verbatim block
      from `afm_zigzag_ribbon.py`'s docstring (now one short paragraph pointing to the
      doc page) and rewrote the doc page's version using the audit's own suggested
      concise text, dropping the review-narration phrasing ("deliberate choice, not an
      oversight", "checked directly", "strictly local effect" — the last also an
      overclaim, since a broadened response need not be strictly local). New version
      states plainly: E=0.05t probes broadened in-gap weight; E=Delta is both the ideal
      edge-state energy and the bulk band edge, so a Gaussian window there also
      overlaps bulk states, making a single-realization comparison less selective.

## Gate 3 — tighten tests, docs, and visuals

- [x] Correct `P_up=(I+sigma_y)/2` -> `(I+sigma_z)/2` in
      `tests/operator_ldos_analytic_stochastic_regression.py`'s comments (the numeric
      expected value was already correct; only the derivation text was wrong).
- [x] Fix `tests/operator_dimension_regression.py`'s module docstring to match its
      actual coverage. Went further than a wording fix: since Gate 1 above added real
      duplicate-index and non-contiguous-index test cases, the docstring's claim is now
      simply true rather than needing to be walked back.
- [x] Isolate `tests/sigma_y_regression.py` and `tests/operator_ldos_regression.py`'s
      outputs in temporary directories, matching the pattern already used in
      `operator_ldos_analytic_stochastic_regression.py`. Made `KITEX`/`KITE_TOOLS`
      absolute paths first (needed once the working directory changes), then wrapped
      each script's build/KITEx/KITE-tools/glob/read sequence in a
      `tempfile.TemporaryDirectory` with `os.chdir`. Both re-verified passing after
      the change; confirmed no stray output files land in the repo anymore.
- [x] Fix the LDOS docstring's `"sz"` example label (`src/kite/__init__.py:690-699`) —
      the implementation requires labels to start with `"l"`, so `"sz"` would actually
      be rejected. Changed to `"lsz"`; also de-capitalized "SAME"/"DIFFERENT" in the
      same paragraph while there.
- [x] Document KITE's actual real-valued random-vector distribution in
      `docs/background/spectral.md`. Verified against `Src/Tools/Random.cpp` first
      rather than trusting the audit's claim: the real instantiation's `initA()` draws
      `(2*dist(rng)-1)*sqrt(3)` with `dist` uniform on `[0,1]` — i.e. continuous uniform
      on `[-sqrt(3),sqrt(3)]` (mean 0, variance 1), not Rademacher; the complex
      instantiation draws `exp(i*2*pi*dist(rng))`, a genuine random phase, matching
      what the page already said for the complex case.
- [x] Fix `afm_zigzag_bands.py`'s "dispersionless" vs. later curvature-explanation
      inconsistency — now says "near-flat edge bands, over their localized momentum
      interval" consistently with the rest of the docstring.
- [x] Fix `afm_zigzag_dos.py`'s capitalization typo ("The Plotted curves" -> "The
      plotted curves" — this had been introduced by an earlier de-capitalization pass
      that missed the mid-sentence case) and colloquial "for real" phrasing (removed).
- [x] Pin `mkdocs-material` in CI instead of always installing latest. Pinned to
      `9.7.6` (the version already confirmed working via the local strict-build check
      earlier in this pass) in both `check-docs` and `deploy-docs` jobs.
- [ ] (Larger, separate effort per audit's own recommendation — not done here) Full
      publication-style plot migration (boxed axes, captions instead of in-figure
      prose, compact E/sigma/averaging-count metadata) — flagged as a follow-up task
      rather than attempted in this pass.

## Not re-opened

- Operator-dimension export-time check, AFM bearded/Klein relabeling, chiral-symmetry
  correction, and CI wiring — all re-confirmed correct/fixed by this audit itself.
