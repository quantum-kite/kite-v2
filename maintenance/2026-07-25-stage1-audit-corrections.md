# Stage 1 audit corrections checklist

Tracks disposition of every finding in `output/pdf/kite_stage1_foundational_audit_2026-07-25.tex`
("KITE Stage 1: Foundational correctness, operator-LDOS, documentation, and change-organization
audit"). Updated as work proceeds; not a substitute for reading the audit itself.

## Decision-summary table (Section 1)

| # | Finding | Status | Notes |
|---|---|---|---|
| 1 | **BLOCKER**: deterministic operator-LDOS contracts with `O^T` (sign bug on Hermitian off-diagonal ops) | ✅ Fixed | `Src/Simulation/Custom/SimulationLDOSOperator.cpp`: `orb_mtx(b,a)` → `orb_mtx(a,b)`. Verified against `sigma_y` (audit's own reproducer) and a full 5-operator regression matrix (I, σz, σx, σy, P↑). See `tests/sigma_y_regression.py`, `tests/operator_ldos_regression.py`. |
| 2 | Spectral background: "numerically exact", missing DOS normalization, non-universal scaling rule | ✅ Fixed | `docs/background/spectral.md`: corrected one-particle diagonalization scaling (polynomial, not exponential), added missing `1/(DR)` in the STE estimator, replaced "numerically exact" with an explicit finite-`M`/kernel/stochastic-variance error-budget statement, qualified the `M`-vs-system-size rule. |
| 3 | Altermagnet page infers zero spin moment from traceless Hamiltonian blocks | ✅ Fixed | `docs/documentation/examples/altermagnet_arpes.md`: replaced with the correct BZ-integration proof (using the already-derived `k_x<->k_y` band relation); fixed the "splitting ∝ d(k)" overclaim (leading-order only) and the "cell doubling always required" AFM comparison (now framed via magnetic space groups). |
| 4 | Markov's inequality described as a relative-error guarantee without the actual theorem | ✅ Fixed | `docs/documentation/examples/markov_local_maps.md`: narrowed the claim to what elementary Markov actually proves; states the paper's stronger per-site bound is not reproduced here. Also fixed "flat and featureless" Bloch-momentum wording in the same file. |
| 5 | `custom_two` conventions public for every `p mod 4`, calibration only recorded for `p=1,2` | ✅ Fixed | Derivation was already general (`i^p` period-4 argument); added `tests/custom_two_mod4_regression.py`, which grounds `p=0,3` to the already-validated `p=2` data via linearity of `Z(E)` in `Gamma_mn` (self-consistency to machine precision). Generalized the units paragraph in `rashba_edelstein_graphene_process.py` beyond the one worked example. |
| 6 | README/tutorial/workflow/installation/example-runner instructions disagree | 🟡 Partially done | Fixed: `README.md` prerequisites (Python≥3.9, no required Pybinding, C++17) + added a verified native-`kite.lattice` quickstart; `docs/documentation/workflow.md` no longer claims the interface is Pybinding-based; `docs/documentation/index.md` tutorial's first-calculation snippet now uses native `kite.lattice` (verified to run) instead of `pybinding.repository.graphene`, and its `tools/build/KITE-tools` path typo is fixed. Not yet done: a full single canonical quickstart page consolidating all of the above in one place; `docs/installation.md` itself not yet fully audited beyond the one Homebrew-installer fix (item 9 below). |
| 7 | `local_chern`/`local_chern_map` removed from docs (C++ path crashes) but still callable in Python | ✅ Fixed | User decision: **remove entirely** — not part of kite-v2. Removed the `local_chern`/`local_chern_map` methods, `get_local_chern`/`get_chern_map` properties (the latter had a pre-existing bug: returned `_local_chern` instead of `_local_chern_map`), storage init, and the `LCM`/`STLCM` HDF5 export block from `src/kite/__init__.py`. Verified the file still parses and `config_system()` still works on an unrelated example. Note: the underlying C++ (`Src/Simulation/SimulationLCM.cpp`, `SimulationSLCM.cpp`, `calc_lcm()`/`calc_st_lcm()` dispatch) was left in place but is now permanently unreachable (no HDF5 trigger data can ever be written) rather than physically deleted — flag if you also want that C++ code excised. |
| 8 | Dirty worktree mixes core implementation, tooling, validation, examples, unrelated work, generated files | ⬜ Not started | User has said hold off on committing; no commit-splitting done yet. |
| 9 | `CLAUDE.md` references 3 missing agent files; internal review narration in public prose | ⬜ Not started | |

## Gate A — operator-LDOS (Section 2)

- [x] `orb_mtx(b,a)` → `orb_mtx(a,b)` fix.
- [x] Hermiticity + shape validation at registration (`_validate_ldos_operators`, wired into `ldos()`/`ldos_map()`).
- [x] Five-operator regression matrix (I, σz, σx, σy, P↑) — `tests/operator_ldos_regression.py`.
- [x] `-1/pi` convention stated explicitly in `ldos()`/`ldos_map()` docstrings (previously the ambiguous `Tr[O*Im G]` shorthand).
- [x] Signed-map limitation stated (rho_O not guaranteed non-negative unless O is positive semidefinite).
- [x] Broadening/seed/disorder-realization alignment requirement now documented directly in `ldos()`'s docstring: exact vs. stochastic-map numerical mismatches aren't evidence of a bug unless broadening and disorder setup are first matched.
- [x] Casing "normalization" resolved as a documentation fix, not a rename: added a code comment in `Simulation.hpp` explaining `calc_LDOS_operators`/`calc_ldos_operators` mirrors the *pre-existing* `calc_LDOS`/`calc_ldos` split (not a new inconsistency), and explicitly why `SimulationLDOSOperator.cpp`/`SimulationLDoSMapOperator.cpp` must NOT be renamed to differ only by case — this repo is edited on a case-insensitive filesystem (macOS default) where that would silently collide into one file (hit once already this session).
- [x] Documented the `operators=` label namespace vs. the positional `l0`-`l9` custom-vertex DSL labels in `ldos()`'s docstring — different namespaces, different failure modes (ValueError vs. silent wrong-operator/out-of-bounds).

## Gate B — mathematical documentation (Section 3, 5)

- [x] One-particle diagonalization scaling corrected (`spectral.md`).
- [x] `1/(DR)` inserted in the STE estimator (`spectral.md`).
- [x] "Numerically exact" replaced with explicit error budget (`spectral.md`).
- [x] `M`-vs-system-size scaling claim qualified (`spectral.md`).
- [x] Altermagnet zero-moment proof corrected (`altermagnet_arpes.md`).
- [x] Markov claim narrowed (`markov_local_maps.md`).
- [x] Time units documented (`time_evolution.md`): `hbar=1` convention, `timestep`↔`energy_scale` relation.
- [x] `time_evolution.md`'s "every other calculation is static linear-response" overclaim fixed.
- [x] All four `custom_two` `p mod 4` residues derived/tested; units paragraph generalized.

## Gate C — usability and CI (Section 4, 6)

- [x] `settings.md` missing comma (syntax bug) fixed.
- [x] `documentation/index.md`: removed dead `calculation.md.md` duplicate reference, fixed `tight_binding.md` path, replaced Pybinding-based first-calculation example with a verified native one.
- [x] `documentation/examples/index.md`: fixed invalid `python -m run_all_examples.py`, fixed old repo URL, clarified the runner only exists under `examples/pybinding/`.
- [x] `examples/README.md`: clarified `run_all_examples.py` location/scope; **fixed** (not just documented) the `npoints`/`num_points` `NameError` in `dccond_phosphorene.py`, removed the now-stale warning note.
- [x] `mkdocs.yml` "Addtitional" typo fixed.
- [x] `docs/installation.md`: retired Ruby Homebrew installer replaced with the current bash installer.
- [x] Six files' stale `quantum-kite/kite` → `quantum-kite/kite-v2` GitHub URLs fixed (verified the referenced commit hashes exist in the shared history first).
- [ ] Symlinked `examples/README.md` relative-link dual-context problem (valid in MkDocs, broken on
  GitHub) — root cause: the reference-link definitions at the bottom of the file
  (`examples/README.md:167-183`, e.g. `[lattice-tutorial]: ../tb_model.md`,
  `[graphene-example]: ../examples/graphene.md`) resolve correctly when MkDocs renders the file at
  its symlinked docs path (`docs/documentation/more_examples/additional_examples.md`), but resolve
  to nonexistent repo-root/examples-root paths when GitHub renders the same physical file at
  `examples/README.md`. User decision (2026-07-25): defer this until actual publish time rather
  than risk destabilizing the mkdocs build now; likely fix is switching just these cross-doc
  reference links to absolute `quantum-kite.com` URLs, but not finalized — revisit before publishing.
- [x] `docs/api/kitex.md` stub expansion (executable syntax, exit behavior, HDF5 read/write contract, rescaling, precision/real-complex instantiation table, decomposition, Python/KITEx/KITE-tools version-coupling caveat) — all facts verified directly against `Src/main.cpp`, not guessed.
- [x] `local_chern`/`local_chern_map` removed from `src/kite/__init__.py` (see item 7 above).
- [x] KITE-tools help-text typos fixed (`.KITE-tools`→`./KITE-tools`, "ommited"→"omitted", "soecified"→"specified", `h5_file.h`→`h5_file.h5`) and the top-level mode list completed (was missing `--ARPES`, `--CustomOne`, `--CustomTwo`; also fixed the stale "four main parameters" count). Rebuilt and verified the corrected `--help` output.
- [ ] `mkdocs build --strict` + link checker in CI.
- [ ] KITEx/KITE-tools numerical smoke test + small fast test set (CTest/Pytest) in CI. (Note: the new `tests/*.py` regression scripts from this pass are not yet wired into any CI/CTest target.)

## Other (Section 7)

- [x] Robotic-phrasing cleanup, first pass: removed the clearest process-narration/appeal-to-authority
  instances. Fixed in `docs/api/kite.md` (two "Verified ... not just asserted" headers/sentences
  rewritten as direct statements of the result; one "not a bug" reassurance rewritten as a plain
  statement of the physical effect). Fixed in `docs/documentation/examples/oam_quadrupole_precession.md`
  (removed the internal-review citation "(checked by `cmt-physicist`)" from public prose — this was
  the sharpest instance of the audit's "internal review narration in public prose" complaint).
  Removed the matching internal-agent-name references from `examples/spin_precession_simple.py` and
  `examples/oam_quadrupole_precession.py` docstrings ("reviewed by cmt-physicist" section headers).
  Surveyed the remaining "not just"/"not merely" hits across `docs/**/*.md` (disorder.md,
  altermagnet_arpes.md, orbital_magnetization.md, rashba_edelstein.md, additional_examples.md,
  performance.md) and judged them factual comparisons rather than narrated verification — left as-is.
  Not yet done: a pass over `examples/*.py` docstrings (only the two `cmt-physicist` citations were
  removed there; the broader "verified"/"exactly" density in example scripts, which the audit's raw
  counts include, has not been surveyed).
- [x] `CLAUDE.md` workflow/agent-file mismatch fixed: `careful-executor` and `structure-auditor`
  are user-level subagent definitions (`~/.claude/agents/`), not part of this repo — the doc's
  claimed paths (`.claude/agents/careful-executor.md` etc.) pointed at files that don't exist in
  the repo's own `.claude/`. Reworded both references to say "user-level definition — not checked
  into this repo" instead of citing a nonexistent in-repo path. `cmt-physicist.md` is unaffected
  (it *is* checked into this repo's `.claude/agents/`).
- [x] `.DS_Store` removal / generated-preview-asset triage: added `.DS_Store` and `/tmp/` to
  `.gitignore`; `git rm --cached` the three already-tracked `.DS_Store` files
  (`.DS_Store`, `docs/.DS_Store`, `docs/documentation/.DS_Store`); deleted the three untracked
  ones from disk (`examples/.DS_Store`, `examples/paper/.DS_Store`,
  `examples/paper/Section_4_E_spintronics/.DS_Store`). `tmp/pdfs/` (13MB of rendered audit-PDF
  pages/LaTeX build artifacts) is scratch output, not source — now gitignored rather than removed
  outright (left on disk in case still needed). Checked `examples/plots/*_preview.png`: this is an
  established tracked pattern in this repo (12 preview PNGs already committed there), so the new
  untracked preview PNGs from this session are normal pending additions, not cleanup targets.
- [ ] Dirty-worktree commit splitting into the audit's recommended dependency order (Section 6.3) — deferred until the user is ready to commit.

## New files added during this pass

- `tests/sigma_y_regression.py` — operator-LDOS transpose-bug regression (audit's own reproducer).
- `tests/operator_ldos_regression.py` — full 5-operator regression matrix.
- `tests/custom_two_mod4_regression.py` — mod-4 self-consistency check for `custom_two`.
