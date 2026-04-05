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
Eigen::Array<T, -1, 1> build_exponential(const double t) {
  static constexpr std::array<std::complex<double>, 4> mi_pow = {{ // (-i)^n
      { 1.0,  0.0},
      { 0.0, -1.0},
      {-1.0,  0.0},
      { 0.0,  1.0}
  }};

  const unsigned N_pols = std::max<unsigned>(4, static_cast<unsigned>(std::ceil(1.5 * t)));
  Eigen::Array<T, -1, 1> moments(N_pols);

  for (int n = 0; n < N_pols; ++n) {
    moments(n) = static_cast<double>(n == 0 ? 1 : 2) * mi_pow[n % 4] * std::cyl_bessel_j(n, t);
  }

  return moments;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_evolwave()
{
  debug_message("Entered Simulation::calc_evolwave\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_evol_wave = false;
    try {
      int dummy_variable;
      get_hdf5<
        int>(&dummy_variable, file, (char *)"/Calculation/evol_wave/Time");
      Global.calculate_evol_wave = true;
    } catch (H5::Exception &e) {
      debug_message("evolwave: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_evol_wave = false;
#pragma omp critical
  local_calculate_evol_wave = Global.calculate_evol_wave;
#pragma omp barrier
  if (local_calculate_evol_wave) {
#pragma omp master
    std::cout << "Calculating time evolution of wave packet.\n";
#pragma omp barrier
    double time;
    int n_measures;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<double>(&time, file, (char *)"/Calculation/evol_wave/Time");
      get_hdf5<int>(&n_measures, file, (char *)"/Calculation/evol_wave/Measurements");

      file->close();
      delete file;
    }
    evolwave(time, n_measures);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::evolwave(double t, int measurements)
{
  debug_message("Entered evolwave\n");
  if constexpr (is_tt<std::complex, T>::value) {
    value_type energy_scale;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
      file->close();
      delete file;
    }
#pragma omp barrier
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type dt   = t / measurements;
    const value_type time = t * energy_scale;
    const value_type step = dt * energy_scale;
    const Eigen::Array<T, -1, 1> coefs = build_exponential(step); // TODO:: Add exponential factor if hamiltonian is shifted

    KPM_Vector<T, D> phi(2, *this);
    Eigen::Array<T, -1, 1> ket(r.Sized);
    Eigen::Array<T, -1, 1> bra(r.Sized);
    Eigen::Array<T, -1, -1> results(r.Sized, measurements + 1); // Store initial wavepacket
    results.setZero();

    h.generate_disorder();
    h.generate_twists();
    phi.initiate_phases();
    phi.set_index(0);
    // phi.initiate_wave_packet2(); // TODO: Implement this
    phi.v.col(0) *= std::sqrt(size);
    phi.Exchange_Boundaries();

    results.col(0) = phi.v.col(0);

    for (int i = 0; i < measurements; ++i) {
      ket.setZero();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }

      results.col(i + 1) = ket;
    }

    store_evolwave(results);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_evolwave(const Eigen::Array<T, -1, -1> &results_)
{
  debug_message("Entered store_evolwave\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.evol_wave.resize(r.Sizet, results_.cols());
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
      Global.evol_wave.row(global.index) = results_.row(local.index);
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/evol_wave/States");
    const std::string name(buffer);
    write_hdf5(Global.evol_wave, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left evolwave\n");
}
