# Custom-vertex rank-one/rank-two fixes: missing EnergyShift, a real Eigen bug, and the "raw velocity" sign convention

**Date:** 2026-07-22
**Scope:** `Src/Simulation/SimulationLDoS.cpp`, `Src/Simulation/SimulationSpectral.cpp`,
`Src/Simulation/Custom/SimulationRankOne.cpp`, `Src/Simulation/Custom/SimulationRankTwo.cpp`,
`Src/Simulation/Simulation.hpp`, `src/kite/__init__.py`, `tools/Src/Custom/customone.{hpp,cpp}` (new),
`tools/Src/Tools/parse_input.{hpp,cpp}`, `tools/Src/Tools/calculate.cpp`,
`examples/haldane_orbital_magnetization.py`.

This note documents three separate, independently-verified issues found while building a new
orbital-magnetization example (`calculation.custom_one()`), and the fix implemented for each. Two are
outright bugs (silently wrong results); the third is a missing generalization (present already for the
built-in conductivity machinery, absent from the newer generic custom-vertex machinery).

---

## 1. `ldos_map()`/`spectral_map()` ignored `EnergyShift` (real bug, fixed)

**Symptom:** for an asymmetric `spectrum_range` (i.e. `energy_shift != 0`), requesting a map at a given real
energy `energy_` silently returned the map for a *different* real energy, off by exactly the shift.

**Root cause:** `Src/Simulation/SimulationLDoS.cpp` and `Src/Simulation/SimulationSpectral.cpp` computed
```cpp
const value_type target = energy_ / energy_scale;
```
missing the `- energy_shift` term that every other calculation type (`ldos()`, `custom_two()`,
`custom_singleshot_two()`) already applies in Python before export
(`(energy - config.energy_shift) / config.energy_scale`, `src/kite/__init__.py`). Invisible whenever
`spectrum_range` happens to be symmetric about zero (the shipped examples, `piflux_ldos_map.py` and
`weyl_spectral_map.py`, both are), which is why it went unnoticed.

**Fix:** read `/EnergyShift` in both files; `target = (energy_ - energy_shift) / energy_scale`.

**Verification:** built a piflux vacancy lattice, requested `energy_=0.0` (the Dirac point) with both a
symmetric and an asymmetric `spectrum_range`. Before the fix, the asymmetric case collapsed to a
featureless map (std=0.013, max=0.15) — matching a control run that *deliberately* mistargeted
`energy_=+shift` in a symmetric config (std=0.012, max=0.11) — proving the bug silently retargeted to the
wrong real energy. After the fix, the asymmetric case matches the correct symmetric result
(std=0.076, max=3.43 vs. std=0.086, max=3.85), while the mistargeted-energy control remains distinctly
different. Regression-checked the two shipped examples (both symmetric `spectrum_range`, so this fix is a
no-op for them) — unaffected.

---

## 2. `store_custom_one`: `.adjoint()` on a vector transposes it; `.conjugate()` doesn't (real bug, fixed)

**Symptom:** every Chebyshev moment beyond `n=0` from `calculation.custom_one()` was silently corrupted, in
every Release build of KITEx.

**Root cause:** `Src/Simulation/Custom/SimulationRankOne.cpp`'s `store_custom_one` computed
```cpp
Global.general_gamma.matrix() += 0.5 * (gamma_.matrix() + gamma_.matrix().adjoint());
```
intending to extract `Re(gamma(n))` per moment (see section 3 for why that's the right operation when the
vertex is Hermitian). `gamma_` is a genuine column vector (`Eigen::Matrix<T,-1,1>`); `.adjoint()` conjugates
**and transposes**, returning a row vector — confirmed directly from Eigen's own headers
(`MatrixBase.h`: `AdjointReturnType = CwiseUnaryOp<scalar_conjugate_op, ConstTransposeReturnType>`). The
resulting Nx1-vs-1xN shape mismatch is:
- **not** caught at compile time — `EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE` passes because one dimension on
  each side is templated `Dynamic` (traced the exact boolean logic in `StaticAssert.h`);
- **not** caught at runtime in KITE's actual Release build — the `eigen_assert` in `CwiseBinaryOp.h` that
  would catch it is compiled out under `-DNDEBUG` (confirmed `CMAKE_CXX_FLAGS_RELEASE = -O3 -DNDEBUG` in
  `build/CMakeCache.txt`).

For `n=0` (a 1x1 "vector") the shapes trivially coincide, so it happened to work — which is exactly why this
went unnoticed: `n=0` alone looked fine.

**Fix:** `.adjoint()` → `.conjugate()` (elementwise, shape-preserving, does not transpose).

**Verification:** before the fix, only the `n=0` moment came out purely real; every other moment retained a
spurious nonzero imaginary part. After the fix, all 512 moments come out with exactly zero imaginary part
(when the vertex is genuinely Hermitian — see section 3 for the case where it isn't).

---

## 3. Raw `"v"` is missing a factor of `i` — a generalization that already existed for built-in conductivities, now added to `custom_one`/`custom_two`

This is not a bug in the sense of "wrong code" — it's a design gap: a correction the original built-in
conductivity machinery already implements, that the newer, more flexible custom-vertex machinery never
inherited.

### The underlying fact

KITE's internal `"vx"`/`"vy"` DSL token, as built in `Src/Hamiltonian/HamiltonianRegular.cpp::build_velocity`
and consumed in `Src/Vector/KPM_Vector2D.cpp`, is literally `[H, r]` (hopping × displacement) — **no factor
of `i`**. The textbook, physically Hermitian velocity operator is `v_phys = (i/hbar)[H, r]`; `[H,r]` alone is
provably **anti-Hermitian** for any Hermitian `H, r` (`[H,r]^dagger = -[H,r]`, a completely general fact).
This does not violate quantum mechanics — the *physical* velocity, once you include its `i/hbar` prefactor,
is still exactly Hermitian, as it must be. It just means KITE's raw internal `"v"` building block is one
factor of `i` away from that physical operator, by design (this keeps the C++ core's arithmetic real-valued
in the common case, deferring any needed compensating `i` to whoever assembles a vertex out of `"v"` tokens).

**Consequence:** whether `Tr[T_n(H) * A]` for a vertex `A` built from `"v"` tokens comes out **real** or
**purely imaginary** depends on the **parity** of how many `"v"` tokens are in `A` — an even count of missing
`i` factors combine to a real sign (`i^2 = -1`, `i^4 = +1`, ...), an odd count leaves a genuine, physical
factor of `i` unresolved. Both cases are physically meaningful; the C++ Hermitization step just needs to know
which component (real, or imaginary) to keep instead of discarding as "noise."

### Where this was already handled correctly

`Src/Tools/Gamma2D.cpp` (the **built-in** `conductivity_dc`/`conductivity_optical`/
`conductivity_optical_nonlinear` machinery) already does this:
```cpp
int num_velocities = 0;
for (int i = 0; i < int(indices.size()); i++) num_velocities += indices.at(i).size();
int factor = 1 - (num_velocities % 2) * 2;   // +1 even, -1 odd
...
gamma = gamma * factor;
...
Global.general_gamma.matrix() += (general_gamma.matrix() + factor * general_gamma.matrix().adjoint()) / 2.0;
```
(Here `.adjoint()` is dimensionally safe — `general_gamma` is a genuine square matrix in this code path, not
a vector, so transposing it is meaningful and intended.) This confirms the original authors were fully aware
of the missing-`i` issue and designed around it — for the built-in conductivity paths.

### Where it was missing

`Src/Simulation/Custom/SimulationRankOne.cpp` (`store_custom_one`, backing the generic `custom_one()`) and
`Src/Simulation/Custom/SimulationRankTwo.cpp` (`store_custom_two`, backing `custom_two()`) always did the
unconditional, `factor=+1`-only version — i.e. always assumed an even/Hermitian case, regardless of how many
`"v"` tokens the caller's vertex actually contained. A vertex with an odd velocity count (like the orbital-
magnetization example below) would have its genuine, physical signal — living entirely in the imaginary part
— silently discarded before it was ever written to the output file, with no way to recover it in
post-processing afterward (verified empirically: sending a raw, odd-velocity-count vertex through the
*unpatched* code gave pure noise, ~1e-4–1e-5 magnitude, vs. genuine signal at ~1e-3–6e-3 once fixed).

### The fix

1. **`src/kite/__init__.py`**: `calculation.custom_one()`/`calculation.custom_two()` now auto-count `"v"`-type
   tokens directly from the vertex's operator strings (`_count_velocity_operators`, splitting each
   dot-separated term on `.` and counting tokens starting with `v`) and write the result as
   `NumVelocities` in the HDF5 export. No manual parameter or compensating `i` in the vertex coefficients is
   required from the user — this mirrors `Gamma2D.cpp`'s own `num_velocities`, just auto-detected from the
   vertex definition instead of an explicit direction list. A vertex mixing terms with inconsistent velocity
   counts (no single well-defined parity) raises `ValueError` rather than guessing.
2. **`Src/Simulation/Custom/SimulationRankOne.cpp`** / **`SimulationRankTwo.cpp`**: `store_custom_one`/
   `store_custom_two` now read `NumVelocities` and apply the identical `factor = 1 - (n%2)*2` logic
   `Gamma2D.cpp` already used, generalizing what was previously a hardcoded `factor=+1`.
3. **`examples/haldane_orbital_magnetization.py`**: the vertex was simplified from a 4-term, manually
   `i`-inserted form (`[[0.5j,"rx.vy"],[-0.5j,"ry.vx"],[-0.5j,"vx.ry"],[0.5j,"vy.rx"]]`) to the raw 2-term
   commutator form with plain real coefficients (`[[1.0,"rx.vy"],[-1.0,"ry.vx"]]`) — `NumVelocities=1` is now
   auto-detected, and the (physically genuine) imaginary-valued signal is correctly preserved rather than
   requiring the caller to manually pre-multiply by `i`.
4. **New KITE-tools mode, `--CustomOne`** (`tools/Src/Custom/customone.{hpp,cpp}`, wired into
   `tools/Src/Tools/parse_input.{hpp,cpp}` and `calculate.cpp`): reconstructs the `custom_one()` spectral
   density `Tr[A * delta(E-H)]` from the raw moments using the exact same Jackson-kernel machinery as
   `--DOS` (`tools/Src/Spectral/dos.cpp`). The one addition: before running the standard real-valued kernel
   sum, it rotates the stored `Gamma` by `-i` when `NumVelocities` is odd (`Gamma` in that case is exactly
   `i*Im(gamma)`, per `store_custom_one`'s corrected Hermitization — multiplying by `-i` recovers the real
   number `Im(gamma)` that the rest of the reconstruction expects). This is the **only** place the `i`
   bookkeeping appears in post-processing; everything else is identical to `--DOS`.

**Why the `i` correction lives here and not in the vertex:** Chebyshev moments themselves don't care whether
the traced operator is Hermitian or anti-Hermitian — `T_n(H)` is computed identically either way. The
correction only matters at the point where the code decides which component (real vs. imaginary) of a
noisy, stochastically-estimated complex number is the genuine signal. That decision belongs in the
Hermitization/reconstruction step (`store_custom_one`/`store_custom_two`, and now `--CustomOne`), not baked
into the vertex's coefficients — keeping the vertex itself a plain, physically direct expression of the
operator being traced.

### Verification performed

- Odd case (`NumVelocities=1`, the orbital-magnetization vertex): confirmed real part is exactly `0.0` and
  the full genuine signal (magnitude ~5e-3, matching the previous manually-`i`-inserted version) lives in the
  imaginary part, both before and after being read back out by the new `--CustomOne` KITE-tools mode.
- Even case (`NumVelocities=2`, Kane-Mele spin-Hall `custom_two()`): confirmed `Γ_mn` is exactly Hermitian
  (`Γ = Γ†` to machine precision) — i.e. completely unchanged behavior from before this change, as expected
  (the even branch of the new `factor` logic reduces to exactly the old, hardcoded behavior).
- Full pipeline (`haldane_orbital_magnetization.py` → KITEx → `KITE-tools --CustomOne`): ran end-to-end,
  produced a finite, sane, non-degenerate spectral density with no `NaN`s.
- Regression: plain `--DOS` on an unrelated lattice (`dos_graphene.py`) still works identically after the
  `parse_input.cpp`/`calculate.cpp` changes needed to wire in `--CustomOne`.

### What this does *not* resolve

This work fixes the **mechanics** of getting a correctly-signed real number out of `custom_one()`/
`custom_two()` regardless of vertex parity. It does **not** by itself settle the **overall physical sign**
of the orbital-magnetization result relative to the reference paper (Vidarte et al., Phys. Rev. B,
doi 10.1103/bhg8-c3rw) — that investigation is still open (see the conversation history / a separate
follow-up note for the k-space cross-check and Chern-number-sign debugging). The fixes here make the
machinery trustworthy enough to *have* that conversation on solid ground; they are not themselves a claim
about the final physical sign of the published example.
