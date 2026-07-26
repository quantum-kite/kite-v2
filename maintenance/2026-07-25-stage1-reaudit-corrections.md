# Stage-1 re-audit corrections checklist

Tracks disposition of every finding in the independent re-audit
`kite_v2_stage1_reaudit_2026-07-25.tex`/`.pdf` (audited commit `795fdf8`). This is a
second-pass review of the corrections made for the first Stage-1 audit
(`maintenance/2026-07-25-stage1-audit-corrections.md`). Updated as work proceeds; not a
substitute for reading the re-audit report itself.

Two of the re-audit's findings were independently re-derived from source (not taken on
faith) before starting work: the missing operator-dimension check, and the bearded/Klein
(not ordinary zigzag) edge geometry of the AFM ribbon lattice — both confirmed correct.

## Executive findings (re-audit's own table)

| # | Finding | Severity | Status | Area |
|---|---|---|---|---|
| 1 | Transpose/sign error in deterministic operator-LDOS is corrected | P1 | ✅ Already fixed (prior pass) | Core math |
| 2 | Operator validation doesn't require matrix dimension == lattice orbital count; deterministic kernel can read outside the matrix | P1 | ✅ Fixed | Core safety |
| 3 | Plain and operator maps run separate propagations and separately regenerate disorder, contradicting "same-run" doc language | P1 | 🟡 Doc fixed; C++ redesign deferred | Core semantics |
| 4 | The ribbon has coordination-one boundary atoms — a bearded/Klein termination, not an ordinary zigzag ribbon | P1 | ✅ Fixed (relabeled) | AFM physics |
| 5 | The stated sublattice chiral symmetry is broken by the staggered mass; finite-width bands are not symmetry-pinned exactly to E=±Δ | P1 | ✅ Fixed | AFM physics |
| 6 | AFM map evaluated off-resonance at E=0.05t, sigma=0.1t, while discussed edge levels lie near E=±0.3t | P2 | ✅ Fixed (documented; on-resonance tried and rejected, see Gate 2) | Example design |
| 7 | Five-operator regression doesn't assert stochastic results for 4 of 5 operators; not wired into CI | P2 | ✅ Fixed | Testing |
| 8 | `custom_two` mod-4 test checks algebraic self-consistency, not end-to-end p=0/p=3 KITE executions | P2 | ⬜ Not started | Testing |
| 9 | Several corrected documentation sections retain contradictions or incorrect API/performance statements | P2 | ✅ Fixed | Documentation |
| 10 | Commit history improved but two commits remain too broad; maintenance checklist is stale | P2 | 🟡 Checklist fixed; broad commits left as-is (see Gate 3) | Organization |
| 11 | OAM/quadrupole-precession example and its nav/assets removed cleanly, no dangling references | P3 | ✅ Confirmed clean | Orbital |

## Gate 1 — operator-LDOS correctness

- [x] Enforce full operator dimension (`shape == (N_orb, N_orb)`) at export time, not just square+Hermitian.
      `config_system()` in `src/kite/__init__.py` now computes `total_orbitals` from the lattice and
      rejects any registered custom operator (any label, used by custom_one/custom_two/ldos/ldos_map/
      gaussian_wave_packet alike) whose shape doesn't match, with a `ValueError` naming the label and
      both shapes. Verified: (1) the existing AFM ribbon example still exports cleanly (no false
      positive), (2) a reproduction of the audit's adversarial case (coupling registered before all
      orbitals were named) is now correctly rejected instead of silently accepted.
- [x] Add negative tests: incomplete, late-grown orbital registrations (`tests/operator_dimension_regression.py`,
      4/4 pass) — covers the two failure modes that actually change matrix *size* (the dimension check's
      concern). Duplicate-index and non-contiguous-index registration are a related but distinct footgun
      (two different names silently aliasing the same physical orbital, or gaps in the index range) that
      the dimension check does not by itself catch, since neither changes `len(_orbital_index_collection)`;
      not covered here, flagged as a separate follow-up rather than claimed as tested.
- [x] Compute plain and operator quantities from the same disorder/stochastic propagation, or explicitly
      redesign/document them as independent calculations with independent cost and randomness.
      Took the documentation path (redesigning the C++ dispatch to share one propagation is a larger,
      separate architectural change, not a quick fix): confirmed via `Src/Simulation/GlobalSimulation.cpp:68-72`
      that `calc_LDOS`/`calc_LDOS_operators` and `calc_ldos`/`calc_ldos_operators` are four independent
      dispatches, and via `SimulationLDoSMapOperator.cpp:146-176` that multiple operator *labels* do share
      one stochastic run (the operator loop is nested inside the vector loop) — so the true picture is:
      cheap to add more labels, not cheap to add operators at all. Documented explicitly rather than fixed
      in code; the shared-propagation redesign itself remains open (noted below).
- [x] Correct the "no extra propagation cost" / "same stochastic run" language in
      `docs/documentation/examples/afm_zigzag_ribbon.md`, and add the same cost note to `ldos()`'s and
      `ldos_map()`'s docstrings in `src/kite/__init__.py`.
- [ ] (Follow-up, not done) Actually redesign the C++ dispatch so plain and operator quantities share one
      disorder/stochastic propagation loop, if the independent-cost/independent-randomness status quo is
      judged unacceptable rather than merely something to document accurately.
- [x] Extend the analytic `H=B*sigma_y` test to all five stochastic operators (closed-form Gaussian-map
      expectations given in the re-audit, Section 3.4) using isolated temporary output directories.
      New `tests/operator_ldos_analytic_stochastic_regression.py`: derives `rho_I`, `rho_sigma_z=0`,
      `rho_sigma_x=0`, `rho_sigma_y`, `rho_P_up=0.5*rho_I` from the two-Gaussian closed form at E=B, runs
      from a `tempfile.TemporaryDirectory` (no more shared-cwd output-name collision risk), and asserts
      each stochastic result against its analytic value within a 6-sigma window of the reported stochastic
      standard error. All 5/5 pass against a real KITEx run.
- [x] Wire the fast operator regressions into CI. Added an "Operator-LDOS regression tests" step to
      `.github/workflows/ci.yml`'s `build-and-test` job (runs on both ubuntu-latest and macos-latest,
      after the existing KITEx smoke test), running `sigma_y_regression.py`,
      `operator_ldos_regression.py`, `operator_ldos_analytic_stochastic_regression.py`, and
      `operator_dimension_regression.py` (all self-contained, ~1-7s each locally). Deliberately did NOT
      wire in `custom_two_mod4_regression.py`: it depends on a `kane_mele_spin_hall-output.h5` fixture
      not generated by any CI step, so it isn't self-contained the same way -- see executive finding #8
      below for the related, larger open item (that script's coverage itself needs strengthening, not
      just its CI wiring).

## Gate 2 — AFM ribbon physics

- [x] Decide: relabel the present geometry as bearded/Klein, or rebuild a true zigzag termination.
      Chose relabeling (rebuilding the lattice is a much larger, riskier change; the edge-localized,
      spin-polarized physics demonstration is equally valid on a bearded/Klein edge). Independently
      re-derived the coordination-1-vs-coordination-2 argument from the actual hopping offsets
      `(0,0),(1,-1),(0,-1)` before trusting the audit's claim. Updated prose in `afm_zigzag_bands.py`,
      `afm_zigzag_ribbon.py`, `afm_zigzag_dos.py`, both `_process.py` scripts, the doc page, `mkdocs.yml`'s
      nav entry, and `examples/README.md`'s two rows. Deliberately did NOT rename any of the
      `afm_zigzag_*.py` filenames themselves (lower-risk: avoids churning cross-references, git history,
      and the already-generated output-file naming scheme) — the geometry is now correctly described in
      every human-readable location without a disruptive rename.
- [x] Remove the "chiral-symmetry-protected"/"pinned exactly"/"100% localized" claims; state the correct
      finite-width edge-energy relation, exact pinning only in the well-localized/infinite-width limit.
      Independently verified numerically (not just taken from the audit) before rewriting: at `Ly=24`,
      the edge band matches Delta to machine precision at k=0, stays within ~1e-7 through most of the
      zone, and only visibly departs approaching where it merges with the bulk continuum near the zone
      boundary. Rewrote `afm_zigzag_bands.py`'s docstring (and its plot title/labels) and the doc page's
      model section accordingly; explicitly noted `{sigma_z, H} != 0` once Delta != 0 (the mass term
      commutes rather than anticommutes with the sublattice operator), so the massless-limit argument is
      not full chiral-symmetry protection of the massive Hamiltonian.
- [x] Replace "a single Dirac cone gaps out" with the two-valley (K, K') description
      (`afm_zigzag_bands.py`'s docstring).
- [x] Evaluate the principal map at E=Delta, or explicitly state that E=0.05t probes broadened in-gap
      tails. Tried on-resonance (E=Delta) first rather than just adding a caveat: the clean-case signal
      is genuinely cleaner there, but a single 5%-vacancy disorder realization's response stopped being a
      small, edge-localized disturbance and grew roughly an order of magnitude through the ribbon
      interior instead (checked directly against the row-averaged data, not assumed) — a real effect from
      part of the zone having additional bands close to Delta that vacancies scatter into more easily,
      not a bug, but a much less legible demonstration for a single disorder draw. Reverted to the
      original off-resonance defaults and instead documented, explicitly and with the numbers behind it,
      why that choice was kept (script docstring, doc page, both figure regions).
- [x] Call the ldos_map output an energy-resolved spin-weighted spectral density, not an "equilibrium spin
      density" or "spin density texture". Fixed in `afm_zigzag_ribbon.py`'s docstring and the doc page's
      figure caption.
- [x] Report stochastic uncertainty when claiming clean-case spin-up/spin-down DOS equality (curves agree
      within stochastic error, not bit-for-bit, in the plotted stochastic data). Fixed in
      `afm_zigzag_dos.py`'s docstring, the in-figure caption text (`afm_zigzag_dos_process.py`), and the
      doc page's figure caption; figures regenerated.
- [x] Add actual literature references for the chiral-symmetry/edge-state claims, or remove the implied
      promise that specific citations exist. Checked first: neither `afm_zigzag_bands.py` nor
      `afm_zigzag_ribbon.py` actually contains any citation (confirmed via grep for arXiv/journal-style
      references) — removed the doc page's dangling "see the script's own docstring for ... literature
      references" promise rather than fabricate citations to match it.

## Gate 3 — documentation and maintenance

- [x] Fix `add_orbital_coupling(row,col,value,label)` -> correct signature `(start,last,c,label)`,
      storing `O[last,start]=c`, in `afm_zigzag_ribbon.md`.
- [x] Attribute the single-`custom_one`-per-file limitation to the Python exporter (serializes only
      `_custom_one[0]`, confirmed at `src/kite/__init__.py:2337-2342` — the same index is read
      regardless of how many times `custom_one()` was called), not to KITE-tools reconstruction.
- [x] Fix workflow.md step 1 still describing/linking a `pb.Lattice` Pybinding build despite the page now
      recommending native `kite.lattice`. Step 1 now leads with `kite.lattice.Lattice`, mentions the
      Pybinding alternative second, and the `[lattice]` reference link now points to the native
      `src/kite/lattice.py` source (matching the fix already applied to `documentation/index.md` in the
      prior audit pass) instead of the Pybinding API docs.
- [x] Fix the new AFM GitHub reference: `/tree/master/examples/...py` -> `/blob/master/examples/...py`
      (file link, not directory link). Found 9 more pre-existing instances of the same tree-vs-blob
      mistake across other example doc pages (not introduced by this AFM work) — out of scope for this
      audit pass, spawned as a separate follow-up task rather than silently expanding scope here.
- [x] Fix real-vs-complex white-noise correlator statement in `docs/background/spectral.md:179-181`
      (can't have both `<chi_i chi_j>=0` and `<chi_i* chi_j>=delta_ij` for real chi). Split into the
      complex-random-phase case (where the two conditions are both correct and distinct, since chi* != chi)
      and the real Rademacher case (where they collapse into a single `<chi_i chi_j>=delta_ij` condition).
- [x] Fix residual Markov/exactness contradiction in `markov_local_maps.md:189-194` (still said "the Markov
      bound" guarantees controlled per-site error, contradicting the page's own earlier, correct explanation
      a few paragraphs up; also replaced "exactness" of the deterministic method with "determinism", since
      it's finite-order, not exact).
- [x] Replace remaining internal-review/process-narration prose in new example docstrings. Fixed the two
      specific quotes the re-audit found ("reviewed and corrected convention" in
      `afm_zigzag_ribbon_process.py`; "confirmed against this KITE build" in `afm_zigzag_dos.py`, which
      also mis-attributed the single-custom_one-vertex limit to KITE-tools instead of the Python exporter
      — fixed alongside), plus "mandatory here, not a default: the whole point of the comparison" (same
      file). Then swept all five `afm_zigzag_*.py`/`_process.py` files for every remaining shouty
      all-caps emphasis word (not just NOT/SAME/EXACTLY): found and fixed ~30 instances (BOTH, LOCAL,
      DEFAULT, BELOW, CLEAN, SPECTRAL, DENSITY, PER, SUBLATTICE, GLOBAL, PROJECTORS, PROPERTY, PLOTTED,
      WHOLE, ALSO, WITHIN, and more), leaving only legitimate acronyms (AFM, KITE, DOS, KPM, DELTA as a
      variable name).
- [x] Update the stale `maintenance/2026-07-25-stage1-audit-corrections.md`: marked items 8 (dirty-tree
      splitting) and 9 (CLAUDE.md agent-path repair) as done (both were completed in later commits after
      that file was last touched), and added a note flagging its `oam_quadrupole_precession` references
      as historical rather than current-state claims, rather than editing them out.
- [ ] (Deliberately not done) Split `579dc36` and `e2c6f23` into narrower commits. Both are already pushed
      to the shared `quantum-kite/master` — splitting them now means rewriting public history (rebase +
      force-push), which this session will not do without the user explicitly requesting it (standing
      git-safety rule, not an oversight).
- [x] Add `mkdocs build --strict` and a link checker to CI. Added a `check-docs` job to
      `.github/workflows/ci.yml` running `mkdocs build --strict` on every push and PR (the existing
      `deploy-docs` job only runs on push to master, so PRs previously got no docs validation at all).
      Ran it locally first: it found one real broken link (`docs/background/index.md`'s `[about]` reference
      pointed at the `../about` directory rather than `../about/index.md`), fixed that, then confirmed a
      clean `--strict` build with zero warnings before wiring it into CI (so this doesn't land as an
      immediately-failing check).
- [ ] Note (no action without a separate decision): the symlinked `examples/README.md` GitHub/MkDocs
      relative-link conflict remains open per prior explicit user decision to defer until publish time.

## Not re-opened

- The core transpose/sign fix (`orb_mtx(a,b)` vs `(b,a)`) — re-confirmed correct by the re-audit itself.
- OAM/quadrupole-precession removal — re-confirmed clean by the re-audit itself.
