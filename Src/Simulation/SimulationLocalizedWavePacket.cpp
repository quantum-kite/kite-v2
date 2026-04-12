#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include <cmath>
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "Loop.hpp"

template <typename T>
T jackson(const int n_, const int polynomials_)
{
  const T arg = M_PI / (polynomials_ + 1);
  const T term1 = (polynomials_ - n_ + 1) * std::cos(arg * n_);
  const T term2 = std::sin(arg * n_) / std::tan(arg);
  return (term1 + term2) / (polynomials_ + 1);
}

template <typename T>
Eigen::Array<T, -1, 1> build_window_here(const T min, const T max)
{
  const T width = max - min;
  const unsigned number_polynomials = std::ceil(256 / width);
  Eigen::Array<T, -1, 1> coefs(number_polynomials);
  coefs(0) =
    jackson<T>(0, number_polynomials) * (std::asin(max) - std::asin(min));
  for (unsigned n = 1; n < number_polynomials; ++n)
    coefs(n) = 2 * jackson<T>(n, number_polynomials) *
               (std::sin(n * std::acos(min)) - std::sin(n * std::acos(max))) /
               n;
  coefs /= M_PI;
  return coefs;
}

template <typename T>
Eigen::Array<std::complex<T>, -1, 1> build_exponential(const T t)
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

template <typename T, unsigned D>
void Simulation<T, D>::calc_localized_wavepacket()
{
  debug_message("Entered Simulation::calc_localized_wavepacket\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_localized_wavepacket = false;
    try {
      int dummy_variable;
      get_hdf5<int>(
        &dummy_variable, file,
        (char *)"/Calculation/localized_wave_packet/Measurements"
      );
      Global.calculate_localized_wavepacket = true;
    } catch (H5::Exception &e) {
      debug_message("localized_wavepacket: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_wavepacket = false;
#pragma omp critical
  local_calculate_wavepacket = Global.calculate_localized_wavepacket;
#pragma omp barrier
  if (local_calculate_wavepacket) {
#pragma omp master
    std::cout << "Calculating time evolution of wave packet.\n";
#pragma omp barrier
    value_type time;
    unsigned n_measures;
    std::array<unsigned, D + 1> pos;
    std::array<value_type, 2> energy_window;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(
        &time, file, (char *)"/Calculation/localized_wave_packet/Time"
      );
      get_hdf5<unsigned>(
        &n_measures, file,
        (char *)"/Calculation/localized_wave_packet/Measurements"
      );
      get_hdf5<unsigned>(
        pos.data(), file,
        (char *)"/Calculation/localized_wave_packet/InitialPos"
      );
      get_hdf5<value_type>(
        energy_window.data(), file,
        (char *)"/Calculation/localized_wave_packet/EnergyWindow"
      );

      file->close();
      delete file;
    }
    localized_wavepacket(time, n_measures, pos, energy_window);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::localized_wavepacket(
  const value_type t,
  const unsigned measurements,
  const std::array<unsigned, D + 1> &pos_,
  const std::array<value_type, 2> &energy_window
)
{
  if constexpr (is_tt<std::complex, T>::value) {
    debug_message("Entered localized_wavepacket\n");
    value_type energy_scale;
    value_type energy_shift;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
      get_hdf5<value_type>(&energy_shift, file, (char *)"/EnergyShift");
      file->close();
      delete file;
    }
#pragma omp barrier
    Coordinates<std::ptrdiff_t, D + 1> global(r.Lt);
    const value_type measure_tau = energy_scale * t / measurements;
    const unsigned divisions = std::ceil(measure_tau / 32);
    const value_type tau = measure_tau / divisions;

    Eigen::Array<T, -1, 1> coefs = build_exponential<value_type>(tau);
    const T arg{0.0, -energy_shift * (tau / energy_scale)};
    coefs *= std::exp(arg);

    KPM_Vector<T, D> phi(2, *this);
    Eigen::Array<T, -1, 1> ket(r.Sized);
    Eigen::Array<T, -1, -1> results(r.Sized, measurements + 1);
    results.setZero();

    h.generate_disorder();
    h.generate_twists();

    phi.initiate_phases();
    phi.set_index(0);
    if constexpr (D == 2)
      global.set({pos_[0], pos_[1], pos_[2]});
    else if constexpr (D == 3)
      global.set({pos_[0], pos_[1], pos_[2], pos_[3]});
    phi.build_site(global.index);
    phi.Exchange_Boundaries();

    if (!(energy_window[0] == 0. && energy_window[1] == 0.)) {
      value_type Emin = (energy_window[0] - energy_shift) / energy_scale;
      value_type Emax = (energy_window[1] - energy_shift) / energy_scale;
      Emin = std::max<value_type>(-1.0, std::min<value_type>(1.0, Emin));
      Emax = std::max<value_type>(-1.0, std::min<value_type>(1.0, Emax));
      const Eigen::Array<value_type, -1, 1> coefs =
        build_window_here<value_type>(Emin, Emax);

      Eigen::Array<T, -1, 1> filtered(r.Sized);
      filtered.setZero();

      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        filtered += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      phi.v.col(0) = filtered;
      phi.empty_ghosts(0);
#pragma omp barrier
#pragma omp master
      Global.soma = 0;
#pragma omp critical
      Global.soma += phi.v.col(0).squaredNorm();
#pragma omp barrier
      phi.v.col(0) /= std::sqrt(Global.soma);
      phi.set_index(0);
      phi.Exchange_Boundaries();
    }
    results.col(0) = phi.v.col(0);
    for (unsigned i = 0; i < measurements; ++i) {
      for (unsigned j = 0; j < divisions; ++j) {
        ket.setZero();
        for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
          phi.cheb_iteration(n);
          ket += coefs(n) * phi.v.col(phi.get_index()).array();
        }
        phi.v.col(0) = ket;
        phi.set_index(0);
        phi.Exchange_Boundaries();
      }
      results.col(i + 1) = phi.v.col(0);
    }
    store_localized_wavepacket(results);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_localized_wavepacket(
  const Eigen::Array<T, -1, -1> &results_
)
{
  debug_message("Entered store_localized_wavepacket\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.localized_wavepacket.resize(r.Sizet, results_.cols());
#pragma omp barrier
  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);
      Global.localized_wavepacket.row(global.index) = results_.row(local.index);
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/localized_wave_packet/States");
    const std::string name(buffer);
    write_hdf5(Global.localized_wavepacket, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left localized_wavepacket\n");
}

#define INSTANTIATE_REAL(type)                                                 \
  template Eigen::Array<type, -1, 1>                                           \
  build_window_here<type>(const type, const type);                             \
  template type jackson<type>(const int, const int);                           \
  template Eigen::Array<std::complex<type>, -1, 1> build_exponential(          \
    const type                                                                 \
  );

INSTANTIATE_REAL(float)
INSTANTIATE_REAL(double)
INSTANTIATE_REAL(long double)
#undef INSTANTIATE_REAL

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
