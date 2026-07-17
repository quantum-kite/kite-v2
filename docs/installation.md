KITE is written in C++ with code optimized for large systems and optimal multithreading performance.
Lattices are specified from Python and turned into the configuration (HDF5) file that [KITEx] reads.
The default, dependency-free way to do this is KITE's own `#!python kite.lattice.Lattice` class
(see [`#!python kite.py`][kitepython]); [Pybinding][pybinding] is also supported as an optional interface
for users who already have lattices defined that way — install it with `#!bash pip install "quantum-kite[pybinding]"`.

The KITE team endeavours to assist researchers run KITE on UNIX-based systems, such as GNU/Linux and Mac OS X.
Thus, feel free to contact any of our team members if you have any queries (contacts can be found at the bottom of the landing page). 

In what follows, we provide detailed installation instructions and additional tips for both Linux and MAC users. 

## 1. Download KITE

First download the source code from our official repository on GitHub [repository][repository]:

``` bash
git clone https://github.com/quantum-kite/kite-v2.git
```

!!! info
   
    Git's installation process for Mac users is outlined in section 2.2.


## 2. Get dependencies

* [HDF5][hdf5] (version 1.8.13 or newer, built with the C++ and high-level APIs)
* [FFTW3][fftw3] (float and double precision — long double is not required, see [Section 5.2][fftw3_common_issue])
* [CMake][cmake] (version 3.9 or newer)
* [gcc][gcc] (version 4.8.1 or newer)
* [h5py][h5py]
* [Pybinding][pybinding] (optional, only if you want to use pybinding-based lattices instead of KITE's native `#!python kite.lattice.Lattice`)

(See detailed instructions below. If you'd rather skip dependency-hunting entirely, jump to
[Section 2.4 (conda)][conda_section] or [Using Docker][docker_section].)

!!! info "Eigen3 is bundled — no separate install needed"

    [Eigen3][eigen3] used to be a common source of cross-machine build issues (different versions/paths
    across `apt`, Homebrew, MacPorts, or a manually unzipped copy). As of this version, a pruned copy of
    Eigen 3.4.0 ships with KITE under `#!bash third_party/eigen3/`, and the build uses it by default — you
    do not need to install Eigen3 yourself. If you specifically want to use your own system Eigen3 instead
    (e.g. to pick up a newer version), pass `#!bash -DUSE_SYSTEM_EIGEN=ON` to `#!bash cmake` in [Section 3][kitex_kitetools].

The compiler **must** support *C++17* (or newer) features and [*OpenMP*][openmp] parallelization. 


To enable KITE's [Gaussian wavepacket propagation][calculation-gaussian_wave_packet] functionality, compile the source code with a recent gcc version
(gcc 8.0.0 or newer).
To check the gcc version, you can use the following command in the terminal:

``` bash
g++ --version
```

### 2.1 For Ubuntu users

Hierarchical Data Format (*HDF5*) is used to store the inputs/outputs of the program:

``` bash
sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
```

KITE also requires FFTW3 (float and double precision):

``` bash
sudo apt-get install libfftw3-dev
```

Calculations on KITE are configured using a Python script, and CMake is required to build the C++ core:

``` bash
sudo apt-get install cmake
```

Then install KITE's Python interface and its dependencies (`#!bash numpy`, `#!bash scipy`, `#!bash h5py`) from
the repository root:

``` bash
pip install -e .
```

If you also want to use pybinding-based lattices (optional — KITE's own `#!python kite.lattice.Lattice`
needs nothing extra):

``` bash
pip install -e ".[pybinding]"
```

### 2.2 For Mac OS X users

The installation of KITE's dependencies on Apple machines is slightly more evolved. We provide below a recipe that has been tested on some Mac OS X systems, but users are encouraged to contact the KITE team shall they encounter any difficulties.  

The *Xcode* command-line tools from Apple Developer are required.  Install these using the terminal:

``` bash
xcode-select --install
```

KITE requires an open-source software package management system like [Homebrew][homebrew] or [MacPorts][ports]. We provide here step-by-step instructions for Homebrew (pointers for MacPorts users are given below). To install HomeBrew, run the following command in the terminal and follow the subsequent instructions provided by software:

``` bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install an up-to-date C++ compiler via Homebrew:

``` bash
brew install gcc
```

Now **close** the terminal window, and open a **new terminal** window.

!!! info
    
    1.  The default directory for Homebrew is */usr/local/bin/*.
        Correct this path if Homebrew was installed in a different directory
    
    2.  In the following sections, replace **n** with the version of gcc installed by Homebrew as given by `#!bash brew info gcc`.


The hierarchical Data Format (*HDF5*) is used to store the inputs/outputs of the program. Install *HDF5* from source, _whilst enforcing the C++17 standard_, using:

``` bash
HOMEBREW_CC=gcc-n HOMEBREW_CXX=g++-n HOMEBREW_CXXFLAGS="-std=c++17" brew install hdf5 --build-from-source
```

!!! info
    
    [MacPorts][ports] users can use the following command:

    ``` bash
    sudo port -v install hdf5 +gcc-n +cxx +hl configure.ldflags="-stdlib=libstdc++" configure.cxx_stdlib="libstdc++" configure.cxxflags="-std=c++17" 
    ```

Install CMake, Python, and FFTW3:

``` bash
brew install python Cmake git fftw
```

Calculations on KITE are configured using a Python script. KITE's own `#!python kite.lattice.Lattice`
needs nothing beyond the packages below; Pybinding is optional (see the note at the top of
[Section 2][get_dependencies]):

!!! warning

    To install the python requirements, you **must** run the Homebrew-python version.
    You can find the Homebrew-python binary at `#!bash /opt/homebrew/bin/python3`.

``` bash
/usr/local/bin/python3 -m pip install -e .
```

To also install pybinding (optional):

``` bash
/usr/local/bin/python3 -m pip install -e ".[pybinding]"
```

Next, download the source code by the command given in section 1. When configuring with CMake in
[Section 3][kitex_kitetools], pass the Homebrew-installed compiler explicitly instead of editing
*CMakeLists.txt*:

``` bash
cmake -DCMAKE_C_COMPILER=gcc-n -DCMAKE_CXX_COMPILER=g++-n ..
```

where **n** is the gcc version installed by Homebrew (as given by `#!bash brew info gcc`).

### 2.3 Verified MacPorts recipe

!!! success "Tested end-to-end"

    Unlike the general recipe above, the steps below are the exact commands used to build and run KITE
    successfully on a Mac (Apple Silicon, MacPorts 2.11.3) — including generating a config file, running
    [KITEx][kitex], and post-processing with [KITE-tools][kitetools]. If you use MacPorts rather than
    Homebrew, this path is the safer bet.

Install a real GCC (not Apple Clang, which is what the plain `#!bash gcc`/`#!bash g++` commands actually
invoke on macOS, and which lacks the `#!bash omp.h` header needed for OpenMP):

``` bash
sudo port install gcc14
```

Install HDF5 **built with the same compiler** and with the C++ and high-level APIs enabled. This last part
matters twice over: KITE requires HDF5's C++ API (`#!bash +cxx`, plus `#!bash +hl` for the high-level API),
which many default HDF5 builds omit — and building it with `#!bash +gcc14` avoids a `libc++`/`libstdc++`
ABI mismatch between HDF5 and KITE at link time, which otherwise shows up as confusing linker errors rather
than a clear "wrong compiler" message:

``` bash
sudo port install hdf5 +cxx +hl +gcc14
```

Install FFTW3 — MacPorts splits this into separate ports per precision; KITE needs the double and float
precision ports (not `fftw-3-long`, see [Section 5.2][fftw3_common_issue]):

``` bash
sudo port install fftw-3 fftw-3-single
```

Since Eigen3 is bundled with KITE (see [Section 2][get_dependencies]), no separate Eigen install is needed.

Install KITE's Python interface (native `#!python kite.lattice.Lattice` needs nothing else; pybinding is
optional, see [Section 2][get_dependencies]):

``` bash
pip install -e .
```

When configuring with CMake (see [Section 3][kitex_kitetools] below), pass the MacPorts gcc14 compiler
explicitly instead of editing *CMakeLists.txt*:

``` bash
cmake -DCMAKE_C_COMPILER=/opt/local/bin/gcc-mp-14 -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-14 ..
```

`#!bash cmake ..` should report `#!bash Found HDF5` with the C++/HL components, `#!bash Found FFTW3`
(float/double), and `#!bash Found OpenMP`, and `#!bash make` should complete with only benign
deprecation warnings (safe to ignore, as noted in Section 3).

### 2.4 Using conda (recommended for reproducibility) { #conda_section }

!!! note "Linux/macOS only — see [Docker][docker_section] on Windows"

    This is intended to sidestep the dependency-hunting issues in sections 2.1–2.3 above (compiler/HDF5
    ABI mismatches, hunting down FFTW3's three separate precision packages, Eigen version drift), and is
    the recommended starting point if you want to avoid those. It's CI-verified end-to-end on both
    `ubuntu-latest` and `macos-latest`. A native (non-Docker) Windows/MSVC build is not supported — MSVC
    rejects several constructs this codebase relies on (`M_PI`, the `and`/`or`/`not` word-operators, an
    unsigned OpenMP loop index, an `Eigen::Map` overload that only resolves on LP64 platforms) that would
    need real portability work across ~20 files to fix, and that work has been deliberately deprioritized
    in favor of [Docker][docker_section], which runs unmodified on Windows via Docker Desktop's WSL2
    backend. If you hit a problem on Linux/macOS with the sequence below, please report it so it can be
    tracked down.

#### 2.4.1 Installing conda itself

If you don't already have conda, the [Miniforge][miniforge] distribution is recommended — it's
conda-forge-first (matching this project's `#!bash environment.yml`) and much lighter than a full Anaconda
install. How you get it depends on what you already have on macOS.

**If you use Homebrew:**

``` bash
brew install --cask miniforge
```

Then initialize your shell (only needed once):

``` bash
conda init "$(basename "${SHELL}")"
```

and open a new terminal window.

**If you use MacPorts:**

MacPorts does not package conda/miniforge (`#!bash port search miniforge` returns no match) — there's
no MacPorts-specific install path. Use the official installer script instead, which is
package-manager-agnostic and installs into its own prefix (`#!bash ~/miniforge3` by default), entirely
separate from MacPorts' `#!bash /opt/local`:

``` bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Follow the prompts (accept the license, default install location is fine). When it asks about running
`#!bash conda init`, either accept (so `#!bash conda activate` works in every new terminal) or decline
if you'd rather activate manually each session — if you decline, run
`#!bash source ~/miniforge3/bin/activate` once per terminal session before the commands below.

#### 2.4.2 Does this conflict with an existing MacPorts/Homebrew install?

No, in practice — verified directly on a machine with MacPorts' own `#!bash gcc14`, `#!bash hdf5`, and
`#!bash fftw` already installed and active. Two things make conda's copies win without any manual PATH
surgery:

1. conda-forge's `#!bash cxx-compiler`/`#!bash c-compiler` packages (already in this project's
   `#!bash environment.yml`) set the `#!bash CC`/`#!bash CXX` environment variables directly to conda's own
   compiler when you activate the environment — CMake and `#!bash make` pick those up ahead of any
   plain `#!bash gcc`/`#!bash g++` found via `#!bash PATH`, so a MacPorts- or Homebrew-installed compiler
   never gets used by accident.
2. The `#!bash cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..` step below explicitly points CMake's
   `#!bash find_package`/`#!bash find_library` at the active conda environment first, so it finds conda's
   HDF5/FFTW3 even when MacPorts or Homebrew have their own separate copies installed.

Confirmed by configuring KITE with conda active on a machine that also has MacPorts' `#!bash gcc14`+`#!bash hdf5`+`#!bash fftw`
installed: CMake reported `#!bash Found HDF5` (version 2.1.0, from `#!bash .../miniforge3/envs/kite/...`),
`#!bash Found FFTW3`, and a Clang 19 compiler — all from the conda environment, none from MacPorts.

The one thing actually worth double-checking (this *is* a real, easy-to-hit gotcha, unlike the imagined PATH
conflict above): if you declined `#!bash conda init`, you must re-run `#!bash source ~/miniforge3/bin/activate kite`
in **every new terminal window** before building — otherwise `#!bash conda`/`#!bash cmake`/`#!bash pip` silently
fall back to whatever your system or MacPorts/Homebrew already provides, with no obvious error. Verify with:

``` bash
which cmake        # should be .../miniforge3/envs/kite/bin/cmake
echo $CC $CXX      # should both be set (not empty)
```

From the `#!bash kite/` directory:

``` bash
conda env create -f environment.yml
conda activate kite
pip install -e .
```

This installs a matched C/C++ compiler toolchain, HDF5, FFTW3, and Python dependencies via conda, then KITE's
own Python interface via `#!bash pip install -e .` (kept as a separate step rather than bundled into
`#!bash environment.yml` itself, since the latter would resolve `#!bash .` relative to wherever
`#!bash environment.yml` happens to be — fine here, but not inside a staged Docker build). Then configure and
build (see [Section 3][kitex_kitetools]) with:

``` bash
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_BUILD_TYPE=Release ..
make
```

No compiler flags or *CMakeLists.txt* edits needed — CMake automatically finds everything inside the active
conda environment. Pybinding is intentionally not part of `#!bash environment.yml` (see the note at the top
of [Section 2][get_dependencies]); add it yourself with `#!bash pip install pybinding` inside the activated
environment if you need it.

## 3. KITEx & KITE-tools
From within the `#!bash kite/` directory (containing *CMakeLists.txt* — the [KITE Python package][kitepython]
should already be installed via `#!bash pip install -e .`, as described in [Section 2][get_dependencies]),
run the following commands:

``` bash
mkdir build
cd build
cmake ..
make
```

!!! info

    Any warnings appearing during the compilation process can typically be ignored.

If these commands have run successfully, you will now find [KITEx][kitex] and [KITE-tools][kitetools]  in the `#!bash kite/build/` directory, which are now ready to use!


## 4. Test KITE

To generate an input file using [KITE's python-interface][kitepython], try one of our examples in the `#!bash kite/examples/` directory:

``` bash
python dos_graphene.py
```

!!! info "This example needs the optional pybinding extra"

    `#!bash dos_graphene.py` builds its lattice with Pybinding, so it needs `#!bash pip install -e ".[pybinding]"`
    (see [Section 2][get_dependencies]). If you'd rather not install Pybinding, try `#!bash shinada_single.py`
    or `#!bash weyl_lt.py` instead — both use KITE's native, dependency-free `#!python kite.lattice.Lattice`.

It creates a file named *graphene_lattice-output.h5* that is used as an input for [KITEx][kitex]:

``` bash
../build/KITEx graphene_lattice-output.h5
```

This first example calculates the density of states (DOS) of pristine graphene.
To obtain the file with the DOS-data, you need to [post-process][kitetools] the output with  

``` bash
../build/KITE-tools graphene_lattice-output.h5
```

which generates the appropriate data file. For more details refer to the [tutorial][tutorial].
 

## 5. Common issues

### 5.1 Using your own Eigen3 instead of the bundled copy

Eigen3 ships with KITE (`#!bash third_party/eigen3/`), so `#!bash cmake ..` should never fail to find it. If
you need a different Eigen3 version for some reason, pass `#!bash -DUSE_SYSTEM_EIGEN=ON` when configuring
(`#!bash cmake -DUSE_SYSTEM_EIGEN=ON ..`); this falls back to CMake's `#!bash find_package(Eigen3)`, in which
case CMake must be able to locate your own [Eigen3][eigen3] install.

### 5.2 FFTW3 "not found", or only some precisions found

KITE links FFTW3's float and double precision libraries — **not** long double, which is intentionally not
required (see below). Some package managers split even these two into separate packages instead of a
single install — notably MacPorts (`#!bash fftw-3`, `#!bash fftw-3-single`, see
[Section 2.3][macports_recipe]) and Debian/Ubuntu, where a plain `#!bash apt-get install libfftw3-dev` is
normally enough, but double-check if `#!bash cmake ..` reports a missing precision.

**Long double is not supported for FFT/spectral calculations, on purpose.** conda-forge's `#!bash fftw`
package (used by [Section 2.4][conda_section]) only ships float and double by default — long double isn't
a default build there, and it's a niche precision requirement in practice — so KITE doesn't require it.
Long-double precision is still fully supported everywhere else in KITE (Hamiltonian, general KPM, etc.);
this is specific to the FFT-based spectral feature. If you genuinely need long-double FFT precision, see
["FFTW3 precision limits"][fftw3_precision_limits] in the code structure docs for what rebuilding that
combination yourself would involve.

### 5.3 `pybinding` fails to build from source

`pybinding` (both on PyPI and conda-forge) is frozen at version 0.9.5 and only ships as a source
distribution for some platform/Python combinations — building its C++ extension can fail depending on your
compiler setup. Since `pybinding` is an optional extra (see the note at the top of
[Section 2][get_dependencies]), you don't need it at all to use KITE — the native
`#!python kite.lattice.Lattice` interface (used by `#!bash shinada_single.py`, `#!bash weyl_lt.py`, and
others in `#!bash examples/`) needs no compiled dependency and always installs cleanly.

## 6. Using Docker

If you'd rather not install any dependencies at all, a `#!bash Dockerfile` is provided that builds a
complete, ready-to-use KITE environment (compiler toolchain, HDF5, FFTW3, and the `#!python kite` Python
package, all via the same [conda-forge environment][conda_section] used in Section 2.4):

``` bash
docker build -t kite .
```

Then run it, mounting your working directory so files you create persist outside the container:

``` bash
docker run -it -v $(pwd):/home/kite/work kite
```

This drops you into a shell inside the container with `#!bash KITEx`, `#!bash KITE-tools`, and
`#!python import kite` all ready to use — no further setup needed. Pybinding is not included in the image by
default (see [Section 2][get_dependencies]); install it yourself inside a running container with
`#!bash pip install pybinding` if you need it.

[repository]: https://github.com/quantum-kite/kite-v2
[eigen3]: https://eigen.tuxfamily.org/
[cmake]: https://cmake.org/
[gcc]: https://gcc.gnu.org/
[h5py]: https://www.h5py.org/
[fftw3]: https://www.fftw.org/
[calculation-gaussian_wave_packet]: api/kite.md#calculation-gaussian_wave_packet
[hdf5]: https://github.com/HDFGroup/
[openmp]: https://gcc.gnu.org/onlinedocs/libgomp/
[homebrew]: https://brew.sh/
[miniforge]: https://github.com/conda-forge/miniforge
[ports]: https://www.macports.org 
[pybinding]: https://docs.pybinding.site/en/stable/install/quick.html
[tutorial]: documentation/index.md
[kitepython]: api/kite.md
[kitex]: api/kitex.md
[kitetools]: api/kite-tools.md
[kitex_kitetools]: #3-kitex-kite-tools
[get_dependencies]: #2-get-dependencies
[conda_section]: #conda_section
[docker_section]: #6-using-docker
[macports_recipe]: #23-verified-macports-recipe
[fftw3_common_issue]: #52-fftw3-not-found-or-only-some-precisions-found
[fftw3_precision_limits]: documentation/code_structure.md#fftw3-precision-limits


