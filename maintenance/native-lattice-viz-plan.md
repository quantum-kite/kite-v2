# Plan: Native lattice/reciprocal-space/band-structure tools (drop remaining pybinding uses)

## Context

`examples/pybinding/` (see its `README.md`) holds the 9 examples that still need pybinding for more than
basic lattice construction. Three of them — `arpes_bilayer.py`, `arpes_cubic.py`, `arpes_tmd.py` — need
exactly two pybinding capabilities beyond `pb.Lattice` itself: `lattice.reciprocal_vectors()` /
`lattice.brillouin_zone()`, and `pb.results.make_path()` for building a k-path through high-symmetry
points. Reproducing these natively (in `kite.lattice`) is one of the "two or three pybinding functions to
reproduce" flagged as future work in that README.

Separately, the user wants to actually **visualize** things KITE currently has no native way to show at
all: the unit cell with its neighbors/hoppings, the Brillouin zone itself (not just use it to pick k-path
points), and a tight-binding band structure along a k-path — useful for sanity-checking a lattice
definition before running an expensive KPM calculation, independent of pybinding.

Scoped by two research passes (not implementation) from `sci-viz` (API surface, what pybinding actually
does, what's simple vs. genuinely complex) and `cmt-physicist` (formula correctness, sign conventions,
consistency with KITE's own established conventions). Findings below are synthesized from both, with one
overstatement from the `sci-viz` pass corrected (noted inline).

**Nothing here is implemented yet — this is the plan for review.**

---

## Where this lives

- **`src/kite/lattice.py`** (no new dependencies — stays importable in headless/HDF5-generation-only
  contexts): add `Lattice.reciprocal_vectors()` and `Lattice.brillouin_zone()` as pure numpy/scipy methods.
  These are reused by non-plotting code too (k-path construction, `H(k)` needs the same reciprocal-space
  facts), so they belong on the data class itself, not bundled with plotting.
- **New `src/kite/visualize.py`** (imports `matplotlib`, `scipy.spatial`): `make_path()`, `plot_unit_cell()`,
  `plot_brillouin_zone()`, `hamiltonian_k()`, `compute_bands()`, `plot_bands()`.

**Correcting one overstatement from the `sci-viz` pass:** it claimed `docs/documentation/examples/spectral_function.md`
"already documents this as forthcoming" for `kite.lattice`. Checked directly — that doc page (written
earlier this session) describes `arpes_tmd.py`'s workflow **as it existed when that example still used
pybinding**; `lattice.brillouin_zone()` there refers to pybinding's real, already-existing method, not a
promise about a future native API. The naming choice below is justified on its own merits (matching
pybinding keeps that doc page accurate with minimal changes if `arpes_tmd.py` is ever converted, and it's
what the user actually wants — see below) — not because it was already decided somewhere.

**Confirmed with the user:** the actual goal is to mirror pybinding's names/functionality specifically so
`arpes_bilayer.py`/`arpes_cubic.py`/`arpes_tmd.py` can eventually move out of `examples/pybinding/` too,
plus get real lattice/Brillouin-zone visualization that doesn't exist natively today.

---

## 1. Reciprocal lattice vectors

**Formula (dimension-agnostic, no 1D/2D/3D branching):** let $A$ be the $m \times D$ matrix of $m \in
\{1,2,3\}$ real-space primitive vectors as rows ($D \ge m$ — e.g. a 2D lattice defined with 3-component
positions). Requiring $\mathbf b_i \in \mathrm{row\,space}(A)$ (the only physical choice when $D > m$) and
$\mathbf a_i \cdot \mathbf b_j = 2\pi\delta_{ij}$:

$$
B = 2\pi\,(AA^T)^{-1}A
$$

Checked against three cases: reduces to the standard $2\pi A^{-T}$ when $m=D$ (matching the convention
already used in `src/kite/__init__.py`'s k-vector fractional/Cartesian conversion, `k_frac = k_cart @
vectors.T / (2π)`); gives the correct zero-transverse-component 1D result; matches the standard
cross-product slab prescription for $m=2, D=3$ (used as a cross-check in tests, not the implementation).

**Do not port the C++ formula literally.** `tools/Src/Spectral/arpes.cpp`'s `2π · vectors.inverse().transpose()`
is only correct because of an *implicit cancellation*: `vectors` is read via a raw `H5Dread` with no
transpose handling into a column-major Eigen array, silently transposing the row-major HDF5 buffer — so
the "extra" transpose in the C++ formula cancels a storage-order transpose, not a mathematical one. Net
result is correct today, but fragile and non-obvious. **Flagged as a separate, out-of-scope-for-now
follow-up**: worth a one-line comment in `arpes.cpp` explaining this, so a future refactor doesn't
"simplify" it into something wrong. Port the *Python-side* relation above, already validated and shipped.

## 2. Brillouin zone

Exact definition: the Voronoi cell of the reciprocal lattice centered at $\mathbf k = 0$. `scipy.spatial.Voronoi`
(Qhull-based) is dimension-generic — 2D and 3D use the same construction call; contrary to the initial
scoping pass, **Voronoi construction itself is not the 2D/3D complexity split**. The real complexity is
downstream, in extracting a *renderable* shape: 2D gives a convex polygon (trivial, ordered vertex list →
`matplotlib.patches.Polygon`); 3D gives a convex polyhedron, and extracting closed, correctly-wound faces
from Qhull's ridge output for `Poly3DCollection` rendering is real, fiddly work (pybinding itself doesn't
support 3D BZ plotting at all — raises `RuntimeError`).

**Decision: implement the Voronoi-cell extraction generically (works for $m=2,3$) from the start — no
extra cost — but only ship *plotting* for 1D (mark $\pm b_1/2$) and 2D (polygon) initially.** 3D BZ
rendering raises a clear "not supported yet" error rather than being attempted half-heartedly. Use enough
candidate reciprocal-lattice points (convergence-checked, not a hardcoded shell radius) so oblique/skewed
cells are still handled correctly.

## 3. k-path (`make_path`)

Pybinding's version is ~25 lines — linear interpolation between consecutive high-symmetry points, tracking
which output indices correspond to the input points (for tick labels later). A plain function returning
something like `(k_points, tick_indices, tick_labels)` is enough — no need for pybinding's custom
`ndarray` subclass.

Two correctness requirements from the physics review, both real (not cosmetic):
- **Point density per segment must scale with actual reciprocal-space Cartesian distance**, not a fixed
  count or fractional-coordinate distance — otherwise apparent band velocity near crossings looks
  distorted between segments of different physical length.
- **Reject zero-length segments** (accidentally repeated consecutive points) rather than dividing by zero.
- **Document the units convention explicitly and pick one**: high-symmetry points supplied by the user are
  Cartesian reciprocal-space coordinates (matching `reciprocal_vectors()`'s output units), not fractional
  — this is exactly the class of unit-convention bug already hit and fixed once this session (ARPES
  k-vectors); don't silently accept both.

## 4. Unit cell + hoppings visualization

Minimal version: scatter each `Sublattice.position` (colored/labeled by name), then for every
`HoppingTerms`, draw a line from `from_sub`'s position to `to_sub`'s position `+ relative_index @ vectors`,
repeated over a small range of neighboring cells (e.g. `n1, n2 ∈ {-1,0,1}`) so boundary-crossing hoppings
are visible. Support 3D via `mplot3d` from the start — `weyl_lt.py` already defines a 3-vector lattice, and
3D scatter+line plotting has no extra geometric complexity (unlike 3D BZ rendering).

## 5. Tight-binding Bloch Hamiltonian + bands — the part most likely to be silently wrong if rushed

**Consistency anchor**: KITE's existing plane-wave state (already shipped, `docs/documentation/examples/spectral_function.md`)
is $|\mathbf k\rangle = \frac{1}{\sqrt N}\sum_{\mathbf r,\alpha} w_\alpha\, e^{i\mathbf k\cdot(\mathbf
r+\mathbf d_\alpha)}|\mathbf r,\alpha\rangle$ — phase carries the **full atomic position** (cell + sublattice
offset), not just the cell index. Any `H(k)` builder must use the same gauge, or its bands won't line up
with KITE's own $A(\mathbf k, E)$ at the same nominal $\mathbf k$.

**Derived (not assumed) sign convention**: Fourier-transforming a `HoppingFamily` entry
(`relative_index=R`, `from_sub=a`, `to_sub=b`, `energy=t_{ab}(R)`, row=to-orbital/col=from-orbital, matching
the existing h.c.-generation code in `__init__.py`) with that same plane-wave phase gives

$$
H_{ba}(\mathbf k) \mathrel{+}= t_{ab}(\mathbf R)\; e^{-i\,\mathbf k\cdot(\mathbf R + \mathbf d_b - \mathbf d_a)}
$$

**with a minus sign in the exponent** — not the naive plus. This is the exact bug class already found and
fixed once this session (commit `e62c839`, the ARPES k-vector fix): a `+` sign is equivalent to evaluating
at $-\mathbf k$, invisible for centrosymmetric bands but visibly mirror-flipped for anything
valley/Rashba/Weyl-asymmetric — get this wrong and it'll look plausible right up until someone checks a TMD
valley or a Weyl point.

**Hermitian-conjugate handling**: `add_one_hopping`'s docstring is accurate — it really does store only one
direction. But the codebase's own established pattern (`__init__.py`'s existing HDF5-export loop) is to
generate the h.c. term automatically at the point of use, not require the caller to add both directions.
**The `H(k)` builder should do the same** — unconditionally add $H_0(\mathbf k) + H_0(\mathbf k)^\dagger$
itself. Trusting users to remember both directions is a quiet, easy-to-miss source of non-Hermitian
matrices (complex eigenvalues don't always look obviously wrong, especially near band touchings). Keep
onsite energies (`Sublattice.energy`) and intra-cell off-diagonal terms (`relative_index=0`,
`from_sub≠to_sub`) logically separate — only Hermitize the latter, never double the onsite block.

**Runtime defense-in-depth**: after building `H(k)` for any sample k, assert it's numerically Hermitian
(`np.allclose(H, H.conj().T)`) and raise rather than silently diagonalizing a broken matrix if not — cheap,
catches wiring/orientation bugs immediately instead of producing subtly-wrong bands.

**API:**
```python
def hamiltonian_k(lattice, k) -> np.ndarray       # Hermitian (n_orb, n_orb) complex, asserts Hermiticity
def compute_bands(lattice, k_path) -> np.ndarray   # eigh() at every k_path point -> (n_k, n_orb) real
def plot_bands(lattice, k_path, ax=None, ylabel="Energy (eV)") -> Axes
```

---

## Implementation order

1. `Lattice.reciprocal_vectors()` — no dependencies.
2. `Lattice.brillouin_zone()` — depends on 1; generic Voronoi extraction, 1D/2D plotting only.
3. `visualize.make_path()` — independent, typically fed `brillouin_zone()` output.
4. `visualize.plot_unit_cell()` — depends only on raw `Lattice` data; can be built/tested first if desired.
5. `visualize.plot_brillouin_zone()` — depends on 2, optionally overlays a path from 3.
6. `visualize.hamiltonian_k()` — depends on `Lattice` data only; verify against `__init__.py`'s h.c.
   generation logic and the derived minus-sign convention above.
7. `visualize.compute_bands()` / `plot_bands()` — depends on 3 and 6.

**Verification before trusting this for anything unfamiliar**: reproduce a known model exactly — graphene's
Dirac cone at K is the standard, cheap sanity check (right energy, right location, right degeneracy). The
existing `examples/pybinding/arpes_tmd.py` 3-band TMD is a second, already-in-repo cross-check once KPM+ARPES
bands can be compared against the same model's tight-binding bands directly.

## Explicitly out of scope for this round

- 3D Brillouin zone polyhedron rendering (construction logic is fine to build generically; rendering isn't).
- Fixing the fragile Eigen row/column-major cancellation in `tools/Src/Spectral/arpes.cpp` (correct today,
  just fragile — worth a comment, not a rewrite, and unrelated to this feature).
- Converting `arpes_bilayer.py`/`arpes_cubic.py`/`arpes_tmd.py` off pybinding themselves — this plan just
  builds the native tools those examples would need; the actual conversion is separate follow-up work.

---

## Addendum — first slice implemented and verified (2D lattice + Brillouin zone)

Per the user's direction to start simple — 2D-only unit-cell/BZ visualization, `H(k)`/bands explicitly
deferred — implemented and independently verified:

- **`Lattice.reciprocal_vectors()`** and **`Lattice.brillouin_zone()`** added to `src/kite/lattice.py`
  (purely additive — `src/kite/__init__.py` confirmed untouched via `git diff`; no existing method on
  `Lattice`/`Sublattice`/`HoppingFamily`/`HoppingTerms` modified). Implements exactly the formulas derived
  above (`B = 2π(AA^T)^{-1}A`; Voronoi cell of the reciprocal lattice, 1D/2D only, `NotImplementedError` for
  3D, fails closed with `RuntimeError` if the extracted region is unbounded rather than returning a wrong
  shape).
- **New `src/kite/visualize.py`**: `plot_unit_cell()` and `plot_brillouin_zone()`, both 2D-only (raise
  `NotImplementedError` otherwise), stable per-sublattice coloring (not matplotlib's default cycle, which
  repeats after 10), plain function-per-figure style matching the rest of the repo's plotting code.
- **Verified independently** (not just trusting the implementing agent's own report) against
  `examples/dos_graphene.py`'s graphene lattice: duality check `a_i·b_j/2π` gives the identity matrix
  exactly; `|b1|=|b2|=29.4987 nm⁻¹` matches the closed-form graphene value; the Brillouin zone comes out as
  a mathematically exact regular hexagon (6 vertices, identical radius, angles exactly 60° apart); the
  rendered figure was viewed directly and shows the correct honeycomb A/B sublattice structure with 3 bonds
  per atom, and a correctly-oriented hexagonal BZ with `b1`/`b2` arrows.

**A real, separate bug found and fixed while committing this**: the new `visualize.py` file was initially
created at `Src/kite/visualize.py` (capital `S`) instead of `src/kite/visualize.py` (lowercase, matching
its siblings `__init__.py`/`lattice.py`/etc.). This isn't a typo — `git status`/`git add` on this repo
kept resolving the new file to the `Src` casing (the pre-existing top-level C++ source directory) no
matter what case was typed in the command, because this repo's `core.ignorecase=true` (git's automatic,
sensible default on macOS's case-insensitive-but-case-preserving filesystem) makes git's working-tree scan
resolve a genuinely new path through the filesystem's already-established directory entry case, which for
the top-level `Src`/`src` collision is `Src` (that directory has existed since the original C++ source
tree, long before this session's `src/kite/` Python package). Confirmed via matching inode numbers that
`Src/` and `src/` are literally the same physical directory on this machine. This does **not** affect
Linux/CI (case-sensitive filesystems correctly split `src/` and `Src/` into two real, separate top-level
directories on checkout — confirmed by `build-and-test (ubuntu-latest)` passing throughout this session
with the existing `src/kite/*` files) — it's specific to adding *new* files under `src/kite/` on this
case-insensitive dev machine. Fixed by staging the file via git plumbing (`git hash-object -w` +
`git update-index --add --cacheinfo` with the explicit correct-case path) rather than `git add`, which
kept re-resolving to the wrong case. **Worth remembering for any future new file added under `src/kite/`
on this machine**: verify with `git status`/`git ls-files` that it landed at the lowercase `src/kite/...`
path, not `Src/kite/...`, before committing.

---

## Addendum — k-path (`make_path`) and predefined lattice templates (graphene, phosphorene, TMDs)

Second slice, per the user's explicit "next stage" request: a k-path helper for the Brillouin zone (item 3
above) and native, pybinding-free equivalents of `pybinding.repository`'s lattice templates, so example
scripts can eventually stop needing pybinding just to get a lattice definition. Since pybinding is actually
installed locally, this was implemented as a direct **port**, not a re-derivation — a `careful-executor`
agent (correctness) built it against pybinding's real source, and a `sci-viz` agent (presentation) then
rendered and visually inspected the result. No `cmt-physicist` involved this round, since there was nothing
to derive — the goal was matching pybinding's existing structure/parameters exactly, not deriving new physics.

**`make_path(*points, step=0.1, point_labels=None)`** added to `src/kite/visualize.py`, wired into
`plot_brillouin_zone`'s existing `k_path` parameter (backward compatible — a plain `(N,2)` array still
works as before). Implements exactly the two correctness requirements from section 3 above: point density
per segment is set by that segment's actual Cartesian length (not fractional coordinates, not a fixed
count), and zero-length segments raise `ValueError` rather than silently dividing by zero. Units convention
(Cartesian reciprocal-space coordinates, matching `reciprocal_vectors()`'s output, not fractional) is
documented explicitly in the docstring, with an explicit callout of the ARPES k-vector unit bug
(`e62c839`) as the reason this isn't silently dual-accepted.

**`src/kite/repository.py`** (new file): `graphene.monolayer(nearest_neighbors=1|2, ...)`,
`phosphorene.monolayer_4band(num_hoppings=2..5)`, `group6_tmd.monolayer_3band(name, ...)` for all 6 group-6
TMD species (MoS2, WS2, MoSe2, WSe2, MoTe2, WTe2). Deliberately mirrors `pybinding.repository`'s
module/function/constant naming (`graphene.a`, `graphene.monolayer`, etc.), per the user's explicit
direction — this is a port, no `pybinding` import anywhere in the file. TMD's 3-orbital metal site is built
as one sublattice per orbital with scalar hoppings (matching the pattern already used natively in
`examples/pybinding/arpes_tmd.py`'s `tmd_monolayer`), since `kite.lattice.Lattice.add_one_sublattice` only
accepts a scalar onsite energy, unlike pybinding's dense per-site onsite matrix.

**Verification** (independently re-checked, not just trusting each agent's own report):
- Every literal constant (graphene `a`/`a_cc`/`t`/`t_nn`, phosphorene's `a`/`ax`/`ay`/`theta`/`phi`/`t1`–`t5`,
  all 6 TMD species' 9-parameter table) diffed exactly (`0.0` difference) against the real, locally-installed
  `pybinding.repository` source — re-ran this diff independently after the agents finished (not just taking
  their word for it): `repository.graphene.a/a_cc/t/t_nn` all compare `True` against
  `pybinding.repository.graphene.a/a_cc/constants.t/constants.t_nn`.
- Independently smoke-tested every template (`graphene.monolayer()` both NN shells, `phosphorene.monolayer_4band()`,
  all 6 `group6_tmd.monolayer_3band(...)` species) — correct sublattice/hopping-term counts (graphene:
  2 subs/3 or 9 hop terms; phosphorene: 4 subs/20 hop terms; each TMD: 3 subs/27 hop terms), all
  `plot_unit_cell`/`plot_brillouin_zone` calls render without error, error-handling paths
  (`make_path` zero-length segment, `graphene.monolayer(nearest_neighbors=3)`) correctly raise `ValueError`.
- `make_path` end-to-end through graphene's actual Γ→K→M→Γ (K/M computed from `reciprocal_vectors()`/
  `brillouin_zone()`, not hardcoded textbook values): K lands exactly on a BZ hexagon corner, M at an edge
  midpoint (both to machine precision), tick_indices index back to the exact input points.
- `sci-viz` rendered and visually inspected the actual figures (not just reasoned about the code) and found
  two real presentation bugs, both fixed in `visualize.py`: (1) `plot_unit_cell` crashed outright on any 2D
  lattice whose sublattice positions carry a 3rd (out-of-plane) component — phosphorene's buckled sites hit
  this immediately — fixed by zero-padding the in-plane cell-translation offset to match the position's own
  dimensionality before adding (a correct top-down projection, not a misrepresentation, since only x/y are
  ever plotted); (2) more seriously, TMD's 3 orbital-sublattices (`Mo_0`/`Mo_1`/`Mo_2`) all sit at the exact
  same in-plane position, and the original same-size-opaque-marker scatter silently let the last-drawn
  orbital hide the other two — the legend showed 3 entries but the figure only ever showed 1. Fixed by
  drawing coincident-position sublattice groups as nested rings (decreasing size, largest first) so every
  orbital stays visible.
- `src/kite/__init__.py` and `src/kite/lattice.py` confirmed untouched (`git diff` shows zero changes) —
  the standing constraint that other tests run against `kite.py` in parallel elsewhere still holds.

**Deliberately deferred** (flagged by the implementing agent, not oversights): graphene's `monolayer_alt`
(same physics, alternative lattice vectors), `monolayer_4atom`/`bilayer` (both need
`Lattice.add_aliases`, which `kite.lattice.Lattice` doesn't implement — a real API gap, not faked around),
and the 3rd-nearest-neighbor `t_nnn` term (explicitly rejected with `ValueError`, not silently ignored).
`hamiltonian_k`/`compute_bands`/`plot_bands` (section 5 above) remain unbuilt — this round was scoped to
`make_path` + lattice templates only. Converting `arpes_bilayer.py`/`arpes_cubic.py`/`arpes_tmd.py` off
pybinding themselves is still separate follow-up work, now unblocked by this landing.

**The same case-sensitivity gotcha recurred** (`src/kite/repository.py`, a new file under `src/kite/`, was
staged by `git status` as `Src/kite/repository.py` again) — same fix as before (`git hash-object -w` +
`git update-index --add --cacheinfo` with the explicit lowercase path), confirmed via `git status --short --
src/` and `git ls-files --others` resolving correctly before committing.

Committed, not yet pushed (in progress alongside a separate CI thread verifying the macOS Bessel-function fix).

---

## Addendum — the 3 remaining ARPES examples converted off pybinding

Per the user's request, once `make_path`/`repository.py` landed: `examples/pybinding/arpes_bilayer.py`,
`arpes_cubic.py`, `arpes_tmd.py` (the only files still needing pybinding beyond basic `pb.Lattice`
construction, per `examples/pybinding/README.md`) converted and moved to the top level of `examples/`.
`pb.Lattice` → `kite.lattice.Lattice`; `pb.results.make_path(...)` → `kite.visualize.make_path(...)[0]`;
`arpes_tmd.py`'s locally-defined `tmd_monolayer()` replaced with a thin wrapper around
`kite.repository.group6_tmd.monolayer_3band()`.

**Verified, not just reviewed**: diffed HDF5 config output (lattice vectors, orbital positions, hoppings,
disorder arrays, weights) between the original pybinding-based version and the converted native version for
all 3 files — bit-identical (max abs diff 0.0) on everything except the `k_vector` path itself, then
independently re-ran all 3 scripts (fresh interpreter, `pybinding in sys.modules` checked explicitly)
to confirm they generate valid output with zero pybinding import at runtime.

**Two things the implementing agent flagged rather than deciding on its own, both resolved**:

1. **A genuine, real bug in the original `arpes_tmd.py`, caught in the process of matching an already-verified
   shared implementation**: its local parameter table used the lattice constant `a` in **angstrom**
   (matching that file's own docstring header), but `kite.repository.group6_tmd`'s parameter table (already
   verified exact-match against `pybinding.repository.group6_tmd`) uses **nm**, like the rest of `kite`.
   Left as a bare substitution, `dk = 0.1` (meant as ≈0.1 Å⁻¹) would have silently become ≈0.1 nm⁻¹ — a real
   ~10x finer/more-expensive k-path sampling than the script ever intended, not caught by the bit-identical
   lattice diff since it only affects the k-path. Fixed by rescaling `dk` to `1.0`; verified independently
   (not just re-trusting the agent) that `distance/step` for each path segment is identical before and after
   to machine precision (13.131003776759846 in both), so the point count is unaffected by the unit change.
   Also fixed the file's stale "Length in angstrom" docstring header to "Length in nm".
2. **A deliberate, accepted design difference between `kite.visualize.make_path` and `pb.results.make_path`**:
   ours uses `ceil(distance/step)` per segment (guarantees at least the requested density, never silently
   drops a segment to 0 points), pybinding's uses `floor` (can undersample short segments). This explains
   every k-path shape mismatch found during verification (e.g. TMD: 27 points old vs. 29 new) — the
   high-symmetry tick points themselves land at numerically identical coordinates either way; only interior
   sampling differs. Decision: keep `ceil` as-is, since it's arguably the more correct/robust behavior (never
   sparser than requested) and was already a deliberate choice when `make_path` was first built — not a
   regression to fix.

Also fixed a real regression the file moves caused: each of the 3 files' `"complete"` CLI mode does
`import run_all_examples as ra` (a driver script that stayed in `examples/pybinding/`) — fixed via
`sys.path.append("pybinding")` in each moved file, and anchoring `run_all_examples.py`'s `_kitex_dir`/
`_kite_tools_dir` to its own file location (`os.path.dirname(os.path.abspath(__file__))`) instead of the
caller's cwd, verified both the direct-invocation and cross-directory-import cases still resolve correctly.

`examples/pybinding/README.md` updated to drop these 3 from its "still need pybinding" list and note that
`kite.repository.graphene.monolayer(...)` also now exists (verified against pybinding) but hasn't yet been
wired into the 3 remaining `dos_*disorder.py` examples — left as the next, not-yet-requested step, along
with `dos_twisted_bilayer.py`'s fundamentally different workflow.

`src/kite/__init__.py`, `lattice.py`, and `repository.py` confirmed untouched throughout.
