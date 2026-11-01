cmake_minimum_required(VERSION 3.9)
set(CMAKE_CXX_STANDARD 17)

# Respect $CC/$CXX (set by conda's compiler-activation packages, or by the
# user directly) instead of hardcoding a compiler. Without an override, CMake
# falls back to its own platform default (e.g. Apple Clang on macOS), which
# may not support the OpenMP flags KITE needs -- pass
# -DCMAKE_C_COMPILER=... -DCMAKE_CXX_COMPILER=... (or set $CC/$CXX) explicitly
# if that's the case for your toolchain.
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DDEBUG=1)
endif()

project(KITE)


find_package(HDF5 1.14.6 COMPONENTS CXX C HL)

# HDF5 is required, not optional, and so are its C++ and high-level components:
# the configuration file, the moments and the results all travel through it, and
# "H5Cpp.h" is included unconditionally. Reporting a missing dependency as a
# status line and configuring anyway turns one clear line into a wall of C++
# errors much later, with the cause scrolled away.
#
# The minimum matches the floor in environment.yml and the one stated in
# docs/installation.md, and it is the lowest version KITE is actually built
# against: the MacPorts route resolves 1.14.6 and the conda route resolves 2.x,
# and both pass the suite. Older releases may well work; none is exercised
# anywhere, and a floor nothing tests is not a supported range.
if(NOT HDF5_FOUND)
    message(FATAL_ERROR
        "KITE requires HDF5 with the C++ and high-level APIs, and CMake could not find it.\n"
        "  Many distributions ship HDF5 without them; KITE needs all three components.\n"
        "  Install a build that has them and point CMake at it:\n"
        "      -DCMAKE_PREFIX_PATH=<prefix>   (or set $HDF5_ROOT)\n"
        "  or use the conda environment, which pins a matched HDF5 alongside the compiler:\n"
        "      conda env create -f environment.yml && conda activate kite\n"
        "  See docs/installation.md for the supported routes.")
endif()

include_directories(${HDF5_INCLUDE_DIR})
set(HDF5_USE_STATIC_LIBRARIES ON)
MESSAGE(STATUS "Hdf5 Library:  ${HDF5_CXX_HL_LIBRARIES}")
MESSAGE(STATUS "Hdf5Hl Library:  ${HDF5_CXX_LIBRARIES}")
MESSAGE(STATUS "Compiler:  ${CMAKE_CXX_COMPILER}")

option(USE_SYSTEM_EIGEN "Use the system-installed Eigen3 instead of the bundled copy in third_party/eigen3" OFF)
if(USE_SYSTEM_EIGEN)
  find_package(Eigen3)
  if(${EIGEN3_FOUND})
    include_directories(${EIGEN3_INCLUDE_DIR})
  else()
    include_directories(~/include/eigen3)
    MESSAGE(STATUS "Couldn't find Eigen3. Looking in ~/include")
  endif()
else()
  MESSAGE(STATUS "Using bundled Eigen3 from third_party/eigen3 (pass -DUSE_SYSTEM_EIGEN=ON to use the system install instead)")
  include_directories(${PROJECT_SOURCE_DIR}/third_party/eigen3)
endif()

# Only float and double precision are required: FFTW3's long-double variant
# isn't built by default by common package managers (e.g. conda-forge's fftw
# package), and long-double FFT/spectral support is a niche need -- see the
# comment on FFTWTraits<long double> in Src/FFT/TraitsFFTW.hpp, which throws a
# clear runtime error if that specific combination is ever actually invoked,
# rather than requiring libfftw3l for every build. Long-double precision is
# still fully supported everywhere else in KITE; this is specific to the new
# FFT-based spectral feature.
find_library(FFTW3_LIB fftw3)
find_library(FFTW3F_LIB fftw3f)
if(FFTW3_LIB AND FFTW3F_LIB)
  MESSAGE(STATUS "Found FFTW3 (float/double): ${FFTW3_LIB} ${FFTW3F_LIB}")
else()
  MESSAGE(STATUS "FFTW3_LIB (double): ${FFTW3_LIB}")
  MESSAGE(STATUS "FFTW3F_LIB (float): ${FFTW3F_LIB}")
  MESSAGE(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
  MESSAGE(FATAL_ERROR "Couldn't find FFTW3 in float and double precision, needed by Src/FFT/ "
    "(see the FFTW3*_LIB/CMAKE_PREFIX_PATH values just printed above for which one(s) are missing and where it looked). "
    "On MacPorts: sudo port install fftw-3 fftw-3-single. On Debian/Ubuntu: sudo apt-get install libfftw3-dev. "
    "On conda: conda install -c conda-forge fftw.")
endif()

set (CORRECT_CODING_FLAGS "")
find_package(OpenMP)
MESSAGE(STATUS "OPENMP:  ${OPENMP_FOUND}")

# OpenMP is required, not optional. `#include <omp.h>` is unconditional in
# Src/Generic.hpp and in twelve translation units under tools/Src/, so a build
# without it does not fall back to a serial path -- every file fails to
# compile. Configuring successfully and failing at compile time pushes a
# missing dependency into a wall of C++ errors, so stop here instead.
if (NOT OPENMP_FOUND)
    message(FATAL_ERROR
        "KITE requires OpenMP, and this toolchain does not provide it.\n"
        "  <omp.h> is included unconditionally, so there is no serial build to fall back to.\n"
        "  On macOS, Apple Clang ships without OpenMP support. Either install a compiler\n"
        "  that has it and point CMake at it:\n"
        "      -DCMAKE_C_COMPILER=<cc> -DCMAKE_CXX_COMPILER=<c++>   (or set $CC/$CXX)\n"
        "  or use the conda environment, which pins a compiler and llvm-openmp together:\n"
        "      conda env create -f environment.yml && conda activate kite\n"
        "  See docs/installation.md for the supported routes.")
endif()

set (CMAKE_C_FLAGS "${CORRECT_CODING_FLAGS} ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CORRECT_CODING_FLAGS} -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CORRECT_CODING_FLAGS} ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")






FILE(GLOB_RECURSE SRCFILES Src/*.cpp)
MESSAGE(STATUS "Files to compile" ${SRCFILES})
add_executable(KITEx ${SRCFILES})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(compile_wp 1)
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0.0)
    MESSAGE( STATUS "Cannot compile gaussian wavepacket")
    set(compile_wp 0)
endif()
add_definitions(-DCOMPILE_WAVEPACKET=${compile_wp})


target_link_libraries(KITEx ${HDF5_CXX_LIBRARIES} ${FFTW3_LIB} ${FFTW3F_LIB})
target_include_directories(KITEx PRIVATE "Src/"
"Src/Tools/"
"Src/Lattice/"
"Src/Hamiltonian/"
"Src/Vector/"
"Src/Simulation/"
)

install (TARGETS KITEx DESTINATION bin)



set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME})
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
set(MACHINE_NAME ${_host_name})

include_directories(${PROJECT_BINARY_DIR}/tools/Src)
FILE(GLOB_RECURSE SRCFILES tools/Src/*.cpp)
add_executable(KITE-tools ${SRCFILES})
string(TIMESTAMP TODAY "(yymmdd) %Y-%m-%d %H:%M:%S")
configure_file (
  "tools/Src/compiletime_info.h.in"
  "tools/Src/compiletime_info.h"
  )

include_directories(${Src})
target_link_libraries(KITE-tools ${HDF5_CXX_LIBRARIES} )
install (TARGETS KITE-tools DESTINATION bin)
