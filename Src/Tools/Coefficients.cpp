#include "Coefficients.hpp"

template <typename T>
T Coefficients<T>::jackson(const int n_, const int polynomials_) const
{
  const T arg = M_PI / (polynomials_ + 1);
  const T term1 = (polynomials_ - n_ + 1) * std::cos(arg * n_);
  const T term2 = std::sin(arg * n_) / std::tan(arg);
  return (term1 + term2) / (polynomials_ + 1);
}

template <typename T>
T Coefficients<T>::gauss_first(
  const unsigned n_,
  const T mu_,
  const T sigma_
) const
{
  const T numerator = n_ * mu_ * sigma_ * sigma_ *
                      (n_ * n_ * sigma_ * sigma_ / (1 - mu_ * mu_) - 3);
  const T denominator = std::pow(1 - mu_ * mu_, -1.5);
  return numerator * denominator;
}

template <typename T>
T Coefficients<T>::gauss_second(
  const unsigned n_,
  const T mu_,
  const T sigma_
) const
{
  const T term_1 = 7 * mu_ * mu_ - 4;
  const T tmp = 1 - mu_ * mu_;
  const T term_2 = 3 - 6 * n_ * n_ * sigma_ * sigma_ / tmp +
                   std::pow(n_ * sigma_, 4) / (tmp * tmp);
  const T denominator = 1 / (24 * tmp);
  return sigma_ * sigma_ * term_1 * term_2 * denominator;
}

template <typename T>
Eigen::Array<T, -1, 1>
Coefficients<T>::build_gaussian(const T energy_, const T width_) const
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
Eigen::Array<T, -1, 1>
Coefficients<T>::build_window(const T energy_, const T width_) const
{
  const T min = energy_ - 0.5 * width_;
  const T max = energy_ + 0.5 * width_;
  const unsigned number_polynomials = std::ceil(64 / width_);
  Eigen::Array<T, -1, 1> coefs(number_polynomials);
  const T diff = std::asin(max) - std::asin(min);
  coefs(0) = jackson(0, number_polynomials) * diff;
  for (unsigned n = 1; n < number_polynomials; ++n)
    coefs(n) = 2 * jackson(n, number_polynomials) *
               (std::sin(n * std::acos(min)) - std::sin(n * std::acos(max))) /
               n;
  coefs /= M_PI;
  return coefs;
}

#define INSTANTIATE_REAL(type) template struct Coefficients<type>;

INSTANTIATE_REAL(float)
INSTANTIATE_REAL(double)
INSTANTIATE_REAL(long double)
#undef INSTANTIATE_REAL
