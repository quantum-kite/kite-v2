# KITE Documentation Project

**Read `GOAL.md` (repo root) at the start of every session.** It's the north star for
this whole project — every stage below should be checked against it, not just against
internal consistency or correctness.

## What this project is
KITE (quantum-kite) is an open-source real-space tight-binding code for quantum transport,
built around Chebyshev-expanded lattice Green's functions (Kernel Polynomial Method).
The C++ core is KITEx, the post-processing tool is KITE-tools, and there's a Python
front-end for building lattices/configurations. Repo: github.com/quantum-kite/kite.

The site (quantum-kite.com) is built by mkdocs from `docs/` at the repo root — NOT from
the separate `quantum-kite/kitedoc` repo, which is old and abandoned. Don't use kitedoc.

We are also comparing against a private fork with reformulated code (added as a separate
remote/clone — check for it before assuming only the public repo is in scope).

## Known documentation gaps (verified against the repo, not guessed)
- `docs/background/tight_binding.md` — conceptual only, no worked Hamiltonian example.
- `docs/api/kitex.md` — a stub (~11 lines). This is the real API gap; `docs/api/kite.md`
  (the Python side) is already solid and shouldn't need a rewrite.
- Domain decomposition (the `divisions` parameter in `kite.Configuration`) is only ever
  mentioned in one-line asides across background/index.md, documentation/settings.md,
  api/kite.md, documentation/disorder.md — never actually explained. No internals/
  contributor-facing section exists at all.
- Domain decomposition is implemented in `Src/Lattice/LatticeStructure.{hpp,cpp}` and
  `Src/Vector/KPM_Vector2D.cpp` / `KPM_Vector3D.cpp`. The boundary/halo region size is
  determined by the Hamiltonian's hopping-neighbor range. There are also code-level
  optimizations here whose *design reasoning* only the original author knows — this
  can't be fully reconstructed from source alone. See workflow below for how to handle it.
- Some doc pages (`docs/about/*.md`, `docs/documentation/more_examples/*.md`) are
  symlinks to root-level files or example READMEs, not standalone files. Preserve this
  pattern rather than duplicating content — edit the README, not a symlink target as if
  it were independent.
- `docs/background/spectral.md` is already solid (KPM, CPGF, stochastic trace evaluation,
  Kubo-Bastin) — don't rewrite it, only extend if a specific gap is found.

## Proposed doc structure
Getting Started (keep) → Background/Theory (keep spectral.md, strengthen tight_binding.md,
add domain-decomposition theory/"why") → API (keep kite.md, flesh out kitex.md) → new
Internals section (domain decomposition implementation + undocumented tricks) → Examples
(keep + expand) → Resources (keep).

## Workflow: three roles, don't skip steps
1. **Plan mode** (built-in) — discuss and scope before touching files. Use this for
   anything non-trivial, especially decisions about doc structure or framing.
2. **careful-executor** (custom subagent, `.claude/agents/careful-executor.md`) — does
   the actual file work. Verifies against source before writing; tags claims as
   `[MECHANICS]` / `[INFERRED REASONING]` / `[NEEDS AUTHOR CHECK]` for anything
   explaining *why* code works a certain way, not just what it does.
3. **structure-auditor** (custom subagent, `.claude/agents/structure-auditor.md`) — read
   only, runs after a multi-file execution, before committing. Checks the *whole* diff
   for organizational coherence: consistency with existing conventions, orphaned files,
   duplicated content, broken nav/cross-references, scope drift, and whether inferred-
   reasoning flags actually survived into the final output.

For the domain decomposition internals page specifically: the design reasoning behind
the boundary sizing and optimizations is not recoverable from code alone. The agreed
fallback is careful-executor drafts its best reconstruction with inferred reasoning
clearly flagged, and the original author reviews/corrects afterward — don't present a
guess as authoritative.

## Physics + visualization review
Two additional subagents exist for this project: `cmt-physicist` (condensed-matter
physics correctness — Chebyshev/KPM formalism, tight-binding, transport; read-only
reviewer, doesn't edit files) and `sci-viz` (matplotlib/scientific visualization —
figure design, presentation, plotting code; full file/execution access).

They can be invoked together on the same figure or physics-bearing script — one
checking correctness, the other checking presentation — producing two independent
reports rather than one blended opinion. Reach for this pairing proactively when a
figure or a physics-bearing claim is genuinely in play, the same as any other
subagent — the user can also request them explicitly by name at any time.

## Style notes
Direct, precise, non-generic. No filler words ("framework", "leveraging", "comprehensive",
"paving the way", etc.) unless they carry real technical content. Prefer LaTeX for any
equations. Don't invent physics, API behavior, or authorship claims — verify or flag
uncertainty instead.
