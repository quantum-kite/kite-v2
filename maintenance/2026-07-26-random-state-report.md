# Random-state audit (2026-07-26)

An external audit (`kite_v2_random_state_and_seed_report_2026-07-26.pdf`, stored outside
this repo) reviewed KITE's random-number architecture. Two separate things came out of
it; only the first is documented/actioned here so far.

## Documented now: `seed_h`/`seed_v` reproducibility

`kite.Configuration(seed_h=..., seed_v=...)` already existed but was undocumented.
`seed_h` seeds Hamiltonian-side randomness (disorder, vacancies, random boundary
twists); `seed_v` seeds the stochastic probe vectors used in trace-based calculations
(e.g. `ldos_map`). Default `0` for either means "don't care" -- KITE asks the OS for a
real random seed and never records it, so the run cannot be reproduced later. Any
positive integer makes that whole run reproducible, since `config_system()` writes the
value into the HDF5 file's `/Seed0`/`/Seed1`, and `KITEx` reads it back from there.
Documented in `docs/api/kite.md`'s `Configuration` entry and in the constructor's own
docstring; added input validation (must be an int in `[0, 2**32)`) so a bad value fails
at construction time instead of silently truncating in the HDF5 `u4` export.

## Open, not yet actioned: shared random state for plain vs. operator-weighted maps

The audit's actual blocker finding: when `ldos_map()`/`ldos()` are called with
`operators=[...]`, the "plain" map and the "operator-weighted" map are computed from
**independent** disorder and probe-vector draws, not the same physical sample --
`SimulationLDoS.cpp` and `SimulationLDoSMapOperator.cpp` (and the analogous exact-LDOS
pair) each call `generate_disorder()` and consume the vector RNG separately. Fixed
seeds make the *whole job* reproducible but do not make these two sub-calculations
share one realization. Empirically: comparing the plain and identity-operator (`O=I`)
maps on the same job showed a relative RMS difference of 0.60 and spatial correlation
of -0.004, i.e. effectively uncorrelated, when they should be identical for `O=I`.

This requires a C++ change (fuse the plain/operator accumulation into one shared
propagation loop, or a checkpoint/restore of the RNG state as a smaller transitional
fix) and is unrelated to the seed-documentation work above. Not scheduled yet.
