# Changelog

All notable changes to KITE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and KITE uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

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

- Four configuration scripts that are shipped as examples or test inputs now build their
  lattices with KITE's own `Lattice` class instead of a third-party modelling package, so
  running them needs nothing beyond KITE's stated dependencies.

### Documentation

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
