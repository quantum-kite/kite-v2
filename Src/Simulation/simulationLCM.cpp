#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "Coefficients.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_lcm()
{
  debug_message("Entered Simulation::calc_lcm\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_lcm = false;
    try {
      int dummy_variable;
      get_hdf5<
        int>(&dummy_variable, file, (char *)"/Calculation/LCM/NumDisorder");
      Global.calculate_lcm = true;
    } catch (H5::Exception &e) {
      debug_message("LCM: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_lcm = false;
#pragma omp critical
  local_calculate_lcm = Global.calculate_lcm;
#pragma omp barrier
  if (local_calculate_lcm) {
#pragma omp master
    std::cout << "Calculating Local Chern.\n";
#pragma omp barrier
    int number_samples;
    value_type energy;
    value_type beta;
    std::array<int, 2> pos;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<
        int>(&number_samples, file, (char *)"/Calculation/LCM/NumDisorder");
      get_hdf5<value_type>(&energy, file, (char *)"/Calculation/LCM/Miu");
      get_hdf5<value_type>(&beta, file, (char *)"/Calculation/LCM/Beta");
      get_hdf5<int>(pos.data(), file, (char *)"/Calculation/LCM/Pos");
      file->close();
      delete file;
    }
    lcm(number_samples, beta, energy, pos);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::lcm(
  const unsigned number_samples_,
  const value_type beta_,
  const value_type energy_,
  const std::array<int, 2> &pos_
)
{
  debug_message("Entered LCM\n");
  value_type energy_scale;
  const value_type norm = 4 * M_PI / r.rLat.determinant();
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
    file->close();
    delete file;
  }
#pragma omp barrier
  KPM_Vector<T, D> phi(2, *this);
  KPM_Vector<T, D> ket(1, *this);
  KPM_Vector<T, D> bra(1, *this);

  bra.initiate_phases();

  const value_type target = energy_ / energy_scale;
  const value_type beta = beta_ * energy_scale;
  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi<value_type>(beta, target);

  Coordinates<std::ptrdiff_t, D + 1> global(r.Lt);
  constexpr unsigned cutt = 2048;
  const unsigned samples = cutt * (number_samples_ / cutt);
  const unsigned partial = samples / cutt;
  Eigen::Array<T, -1, 1> results(partial * r.Orb);

  for (unsigned disorder = 0; disorder < samples; disorder += partial) {
    results.setZero();
    for (unsigned p = 0; p < partial; ++p) {
      h.generate_disorder();
      h.generate_twists();
      for (unsigned orb = 0; orb < r.Orb; ++orb) {
        phi.set_index(0);
        global.set({r.Lt[0] / 2 + pos_[0], r.Lt[1] / 2 + pos_[1], orb});
        phi.build_site(global.index);
        phi.Exchange_Boundaries();
        ket.v.col(0).setZero();

        for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
          phi.cheb_iteration(n);
          ket.v.col(0) += phi.v.col(phi.get_index()) * coefs(n);
        }
        bra.v.col(0) = ket.v.col(0);
        ket.mult_position(1, &ket);

        phi.set_index(0);
        phi.v.col(0) = ket.v.col(0);
        phi.cheb_iteration(0);
        ket.v.col(0) = 0.5 * phi.v.col(phi.get_index());
        for (unsigned n = 1, N = coefs.size(); n < N; ++n) {
          phi.cheb_iteration(n);
          ket.v.col(0) -= phi.v.col(phi.get_index()) * coefs(n);
        }
        ket.mult_position(0, &ket);
        bra.empty_ghosts(0);
        for (std::size_t ii = 0, II = r.Sized; ii < II; ii += r.Ld[0]) {
          const auto tmp = bra.v.col(0)
                             .segment(ii, r.Ld[0])
                             .dot(ket.v.col(0).segment(ii, r.Ld[0]));
          results(orb + p * r.Orb) += norm * tmp;
        }
      }
    }
    store_lcm(results);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_lcm(const Eigen::Array<T, -1, 1> &results_)
{
  debug_message("Entered store_gamma\n");
#pragma omp master
  {
    Global.general_gamma.resize(results_.size(), 1);
    Global.general_gamma.setZero();
  }
#pragma omp barrier
#pragma omp critical
  Global.general_gamma += results_;
#pragma omp barrier
#pragma omp master
  {
    Global.buffer.insert(
      Global.buffer.end(), Global.general_gamma.data(),
      Global.general_gamma.data() + Global.general_gamma.size()
    );
    Eigen::Array<T, -1, -1> tmp = Eigen::Map<
      Eigen::Array<T, -1, -1>>(Global.buffer.data(), Global.buffer.size(), 1);
    H5::H5File file(name, H5F_ACC_RDWR);

    const std::string name_dataset = "/Calculation/LCM/Marker";
    if (H5Lexists(file.getId(), name_dataset.c_str(), H5P_DEFAULT))
      file.unlink(name_dataset);
    write_hdf5(tmp, &file, name_dataset);
  }
#pragma omp barrier
  debug_message("Left store_gamma\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
