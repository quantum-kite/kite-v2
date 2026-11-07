# Changelog

All notable changes to KITE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and KITE uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `import kite` printed `SyntaxWarning: invalid escape sequence '\\i'` on every import. A
  docstring documenting the `twisted` boundary option used the mathematical `\\in` without an
  `r` prefix. The docstring is now raw; its text is unchanged. Python has said this warning
  becomes a `SyntaxError`, so the package would eventually have failed to import.

- The example sweep's own guard stopped running. `verdict()` in
  `tests/run_all_examples.py` grew a fourth return value when deferred failures were
  introduced, and `tests/sweep_exit_status_regression.py` kept unpacking three, so five
  of its cases raised before asserting anything. The gate was widened in the same change
  that disabled the check on it. The exit decision now lives in one function,
  `exit_status`, which both the sweep and its guard call, so neither restates the rule.

- `examples/ldos_graphene.py` asked for a spectral range of `[-2.8, 2.8]` eV while using
  `t = 2.8` eV. Graphene's bandwidth is ±3t = ±8.4 eV, so the range was too narrow by a
  factor of three, the Chebyshev recursion diverged, and the example wrote local densities
  of order 1e41 while still exiting 0. The range now brackets the spectrum with margin.
- `tools/tests/next7/config.py` could not run at all: `Calculation.arpes` takes a required
  `weight` argument, one entry per orbital, and the script passed none.
- `examples/spin_precession_simple_process.py` failed to read its own output under NumPy 2.
  `timestep` is stored as a length-1 array, and converting one of those to a Python scalar
  has been an error since NumPy 2.
- `tests/run_all_examples.py` crashed while reporting an example that had timed out,
  instead of reporting it: the timeout carries undecoded output buffers even when the
  subprocess call asked for text.
- The build routes that are not CMake — the plain `makefile` and `compile_script.sh` — did
  not build the Gaussian wave-packet propagation, so the same source produced binaries
  with different capabilities depending on how it was compiled.

### Changed

- The three Weyl optical-conductivity examples (`weyl_lt.py`,
  `optcond_t_symmetric_cubic_weyl_sm.py` and the `pybinding` variant) ran for over half an
  hour each at their published size. They now use a lattice and moment count that finish in
  seconds on a laptop, and each carries a comment saying these are reduced parameters, which
  two to raise to see the physics, and in what order — `num_moments` dominates, because the
  cost is quadratic in it.
- Four configuration scripts that are shipped as examples or test inputs now build their
  lattices with KITE's own `Lattice` class instead of a third-party modelling package, so
  running them needs nothing beyond KITE's stated dependencies.

- `tests/run_all_examples.py` deleted each example's working directory without checking
  what `--outdir` pointed at, so a mistyped path would have removed directories under it.
  It now marks its own output root and refuses any directory holding content it did not
  create. It also wrote absolute paths into each `run.log`, putting the name of whoever
  ran it into a file that `--outdir` can place anywhere; the repository and home prefixes
  are now replaced before the log is written.

### Removed

- `examples/cond_sum.py` and `examples/process_conductivities.py`, two Python helpers that
  summed and post-processed conductivity output. `KITE-tools` performs both, so the scripts
  duplicated the shipped path in a second, unverified implementation.

- `KITE-tools --info` printed nothing. It is advertised in `--help` as reporting the
  compilation details, and the parser accepted it, but the function behind it had its
  whole body commented out while the generated header it needs was being built and was
  already on the include path. It now prints the machine, operating system and
  configuration date, in both the long and short forms.

### Documentation

- The header of `examples/ldos_graphene.py` described the example as building a
  near-decoupled resonant impurity with severed bonds, producing a sublattice-polarised
  LDOS pattern. The script does not do that: the bond entries add to the clean hopping
  rather than replacing it, so a value of `-0.0` changes nothing, and the output is
  spatially flat to within 2%. The header now states what the example is — local density
  of states with structural disorder — and nothing more. No behaviour changed.

- The KITE-tools output table said `dos.dat` carries a third column with the imaginary part
  of the density of states, and that single-shot DC conductivity carries one too. Neither
  does: both files have two columns, and the single-shot file is `singleshot_condDC.dat`,
  not a name derived from the HDF5 file. Verified against output the examples produced, not
  only against the source. The other three entries claiming an imaginary part are correct.

- `docs/installation.md` names CMake as the main build route and describes the two
  alternatives, each with what it does and does not cover; it also gains a section on how
  to run the verification suite and a statement of what continuous integration actually
  checks, which is narrower than the suite.
- The minimum GCC version in the dependency list was still 4.8.1, which has not been
  sufficient since the code moved to C++17. It now says 8.0.0 or newer.
- `docs/documentation/settings.md` says what a too-narrow `spectrum_range` actually does:
  both binaries print their usual completion message and exit 0, and the output file is
  filled with `nan`. Nothing in the run reports the failure.
- `docs/api/kite.md` states what a KITEx run is reproducible from — the two seeds **and**
  `divisions`, because each seed is mixed with the OpenMP thread index. Quoting a seed
  alone does not let anyone reproduce a figure.
- `make_pybinding_model` is documented at the import path it actually has.
- The five post-processing keywords that `KITE-tools` accepts but no page mentioned are
  now documented, with the options each one takes.
- The `Calculation` reference points at the method that derives the expansion order for
  the two map functions, which take no `num_moments` of their own.

## [1.1.0]

Earlier releases are not covered by this file.
