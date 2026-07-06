#ifndef GLOBAL_FFT_H_
#define GLOBAL_FFT_H_
#include <complex>

template <typename T>
struct GlobalFFT {
  std::complex<T> *in = nullptr;
  std::complex<T> *out = nullptr;

  GlobalFFT() = default;
  ~GlobalFFT();
  void allocate(const std::size_t total_size);
};

#endif
