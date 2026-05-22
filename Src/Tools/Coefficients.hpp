#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>

namespace Coefficients {
template <typename T>
T jackson(const int n_, const int polynomials_)
{
  const T arg = M_PI / (polynomials_ + 1);
  const T term1 = (polynomials_ - n_ + 1) * std::cos(arg * n_);
  const T term2 = std::sin(arg * n_) / std::tan(arg);
  return (term1 + term2) / (polynomials_ + 1);
}

template <typename T>
T gauss_first(
  const int n_,
  const T mu_,
  const T sigma_
) {
  const T numerator = n_ * mu_ * sigma_ * sigma_ *
                      (n_ * n_ * sigma_ * sigma_ / (1 - mu_ * mu_) - 3);
  const T denominator = std::pow(1 - mu_ * mu_, -1.5);
  return numerator * denominator;
}

template <typename T>
T gauss_second(
  const int n_,
  const T mu_,
  const T sigma_
) {
  const T term_1 = 7 * mu_ * mu_ - 4;
  const T tmp = 1 - mu_ * mu_;
  const T term_2 = 3 - 6 * n_ * n_ * sigma_ * sigma_ / tmp +
                   std::pow(n_ * sigma_, 4) / (tmp * tmp);
  const T denominator = 1 / (24 * tmp);
  return sigma_ * sigma_ * term_1 * term_2 * denominator;
}

template <typename T>
Eigen::Array<T, -1, 1> build_gaussian(const T energy_, const T width_)
{
  const unsigned number_polynomials = std::ceil(6.0 / width_);
  Eigen::Array<T, -1, 1> coefs(number_polynomials);
  coefs(0) = 1 - gauss_second(0, energy_, width_);

  for (int n = 1; n < number_polynomials; ++n) {
    const T gaussian =
      std::exp(-0.5 * n * n * width_ * width_ / (1 - energy_ * energy_));
    const T cossine = std::cos(n * std::acos(energy_));
    const T sine = std::sin(n * std::acos(energy_));
    coefs(n) = 2 * gaussian *
               (cossine * (1 - gauss_second(n, energy_, width_)) -
                0.5 * sine * gauss_first(n, energy_, width_));
  }
  const T prefactor = 1 / (M_PI * std::sqrt(1 - energy_ * energy_));
  coefs *= prefactor;
  return coefs;
}

template <typename T>
Eigen::Array<T, -1, 1> build_window(const T center_, const T width_)
{
  const T min = center_ - 0.5 * width_;
  const T max = center_ + 0.5 * width_;
  const unsigned number_polynomials = std::ceil(64 / width_);
  Eigen::Array<T, -1, 1> coefs(number_polynomials);
  const T diff = std::asin(max) - std::asin(min);
  coefs(0) = jackson<T>(0, number_polynomials) * diff;
  for (unsigned n = 1; n < number_polynomials; ++n)
    coefs(n) = 2 * jackson<T>(n, number_polynomials) *
               (std::sin(n * std::acos(min)) - std::sin(n * std::acos(max))) /
               n;
  coefs /= M_PI;
  return coefs;
}

template <typename T>
Eigen::Array<std::complex<T>, -1, 1> build_cplx_exp(const T t)
{
  static constexpr std::array<std::complex<T>, 4> mi_pow = {
    {{1.0, 0.0}, {0.0, -1.0}, {-1.0, 0.0}, {0.0, 1.0}}
  };
  const unsigned N_pols =
    std::max<unsigned>(32, static_cast<unsigned>(std::ceil(2 * t)));
  Eigen::Array<std::complex<T>, -1, 1> moments(N_pols);

  moments(0) = mi_pow[0] * static_cast<T>(std::cyl_bessel_j(0, t));
  for (unsigned n = 1; n < N_pols; ++n)
    moments(n) = mi_pow[n % 4] * static_cast<T>(2.0 * std::cyl_bessel_j(n, t));
  return moments;
}

template <typename T>
Eigen::Array<T, -1, 1> build_fermi(const T beta_, const T mu_)
{
  const unsigned N = std::ceil(3.0 * beta_);
  Eigen::Array<T, -1, 1> x(N);
  Eigen::Array<T, -1, 1> f(N);
  Eigen::Array<T, -1, 1> ac(N);
  Eigen::Array<T, -1, 1> coef(N);
  for (unsigned i = 0; i < N; i++) {
    x(i) = std::cos(0.5 * M_PI * (2 * i + 1) / N);
    if (beta_ * (x(i) - mu_) > 200)
      f(i) = 0.0;
    else if (beta_ * (x(i) - mu_) < -200)
      f(i) = 1.0;
    else
      f(i) = 1.0 / (1 + std::exp(beta_ * (x(i) - mu_)));
    ac(i) = std::acos(x(i));
  }
  for (unsigned i = 0; i < N; i++) {
    T sum = 0.0;
    for (unsigned j = 0; j < N; j++)
      sum += std::cos(i * ac(j)) * f(j);
    coef(i) = sum * 2.0 / N;
  }
  coef(0) *= 0.5;
  return coef;
}

}

#endif
