#ifndef FFTW_TRAITS_H_
#define FFTW_TRAITS_H_
#include <fftw3.h>
#include <stdexcept>

template <typename T>
struct FFTWTraits;

template <>
struct FFTWTraits<float> {
  using cplx = fftwf_complex;
  using plan = fftwf_plan;

  static void *malloc(std::size_t size) { return fftwf_malloc(size); }
  static void free(void *ptr) { fftwf_free(ptr); }

  template <typename... Args>
  static plan plan_many_dft(Args &&...args)
  {
    return fftwf_plan_many_dft(std::forward<Args>(args)...);
  }

  static void execute(plan p) { fftwf_execute(p); }
  static void destroy_plan(plan p) { fftwf_destroy_plan(p); }
};

template <>
struct FFTWTraits<double> {
  using cplx = fftw_complex;
  using plan = fftw_plan;

  static void *malloc(std::size_t size) { return fftw_malloc(size); }
  static void free(void *ptr) { fftw_free(ptr); }

  template <typename... Args>
  static plan plan_many_dft(Args &&...args)
  {
    return fftw_plan_many_dft(std::forward<Args>(args)...);
  }

  static void execute(plan p) { fftw_execute(p); }
  static void destroy_plan(plan p) { fftw_destroy_plan(p); }
};

// KITE only links the float and double precision FFTW3 libraries (see the
// FFTW3 detection block in the top-level CMakeLists.txt) -- the long-double
// variant is not built by default by common package managers (e.g.
// conda-forge's `fftw` package) and is a niche precision requirement, so it
// isn't required for a standard build. Long-double precision is still fully
// supported everywhere else in KITE (Hamiltonian, general KPM, etc.) -- this
// restriction is specific to the FFT-based spectral feature. Anyone who
// genuinely needs long-double FFT/spectral support must build FFTW3 from
// source with `--enable-long-double`, link libfftw3l themselves, and replace
// the throwing calls below with real fftwl_* calls (mirroring the float/double
// specializations above).
template <>
struct FFTWTraits<long double> {
  using cplx = fftwl_complex;
  using plan = fftwl_plan;

  [[noreturn]] static void not_supported()
  {
    throw std::runtime_error(
      "FFT/spectral calculations with long-double precision are not supported "
      "in this build: KITE only links FFTW3's float and double precision "
      "libraries by default. Rebuild FFTW3 from source with "
      "--enable-long-double, link libfftw3l, and update Src/FFT/TraitsFFTW.hpp "
      "if you need this combination."
    );
  }

  static void *malloc(std::size_t) { not_supported(); }
  static void free(void *) { not_supported(); }

  template <typename... Args>
  static plan plan_many_dft(Args &&...)
  {
    not_supported();
  }

  static void execute(plan) { not_supported(); }
  static void destroy_plan(plan) { not_supported(); }
};

#endif
