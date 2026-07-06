#ifndef FFTW_TRAITS_H_
#define FFTW_TRAITS_H_
#include <fftw3.h>

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

template <>
struct FFTWTraits<long double> {
  using cplx = fftwl_complex;
  using plan = fftwl_plan;

  static void *malloc(std::size_t size) { return fftwl_malloc(size); }
  static void free(void *ptr) { fftwl_free(ptr); }

  template <typename... Args>
  static plan plan_many_dft(Args &&...args)
  {
    return fftwl_plan_many_dft(std::forward<Args>(args)...);
  }

  static void execute(plan p) { fftwl_execute(p); }
  static void destroy_plan(plan p) { fftwl_destroy_plan(p); }
};

#endif
