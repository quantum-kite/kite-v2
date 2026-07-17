#ifndef BESSEL_J_H_
#define BESSEL_J_H_

#include <cmath>
#include <vector>

// Portable replacement for std::cyl_bessel_j (the C++17 "special math
// functions" library). That standard-library feature is implemented by
// libstdc++ (GCC, used on Linux) but NOT by libc++ (Clang, used on macOS --
// both Apple's toolchain and conda-forge's -- and this was confirmed to fail
// to compile there with "no member named 'cyl_bessel_j' in namespace
// 'std'"). These functions were left optional in the C++17 standard for
// exactly this reason. Vendoring a small, portable implementation here
// avoids depending on a standard-library feature with inconsistent platform
// support, matching how Eigen is already vendored in this project.

// Fills out[0..n_max] with the cylindrical Bessel function of the first
// kind, J_k(x), for integer order k = 0, 1, ..., n_max and real argument x.
//
// Uses Miller's algorithm (downward recurrence): the three-term recurrence
// relation J_{k-1}(x) + J_{k+1}(x) = (2k/x) J_k(x) is numerically unstable
// when recurred upward (increasing k) for fixed x, since J_k(x) decays
// rapidly once k exceeds x and upward recurrence amplifies rounding error.
// It is stable when recurred downward from an order well above both n_max
// and x, starting from arbitrary (unnormalized) seed values; the true,
// correctly-scaled values are then recovered via the normalization identity
// J_0(x) + 2 * sum_{k=1}^{infinity} J_{2k}(x) = 1.
// (Standard method; see e.g. Abramowitz & Stegun 9.1.27, or any numerical
// methods reference covering Bessel function recurrence.)
inline void cyl_bessel_j_series(unsigned n_max, double x, double *out)
{
  if (x == 0.0) {
    out[0] = 1.0;
    for (unsigned k = 1; k <= n_max; ++k)
      out[k] = 0.0;
    return;
  }

  const double ax = std::fabs(x);
  const bool negative_x = x < 0.0; // J_k(-x) = (-1)^k J_k(x)

  // Starting order for the downward recurrence: must be comfortably above
  // both n_max and ax so the (arbitrary) seed values have decayed to
  // negligible relative influence by the time the recurrence reaches n_max.
  unsigned m = n_max + static_cast<unsigned>(std::sqrt(40.0 * (n_max + 1)) + ax) + 10;
  if (m % 2 != 0)
    ++m; // even, so the normalization sum below cleanly covers k = 2, 4, ..., m

  std::vector<double> j(m + 2, 0.0);
  j[m] = 1.0; // arbitrary nonzero seed; overall scale fixed by normalization below

  for (unsigned k = m; k >= 1; --k) {
    j[k - 1] = (2.0 * static_cast<double>(k) / ax) * j[k] - j[k + 1];
  }

  double sum = j[0];
  for (unsigned k = 2; k <= m; k += 2)
    sum += 2.0 * j[k];

  for (unsigned k = 0; k <= n_max; ++k) {
    double value = j[k] / sum;
    if (negative_x && (k % 2 == 1))
      value = -value;
    out[k] = value;
  }
}

#endif
