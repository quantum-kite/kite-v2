#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_
#include <Eigen/Dense>

template <typename T>
struct Coefficients {
  T jackson(const int, const int) const;
  Eigen::Array<T, -1, 1> build_window(const T, const T) const;
  T gauss_first(const unsigned, const T, const T) const;
  T gauss_second(const unsigned, const T, const T) const;
  Eigen::Array<T, -1, 1> build_gaussian(const T, const T) const;
};

#endif
