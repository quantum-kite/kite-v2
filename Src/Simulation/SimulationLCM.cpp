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
#include "Loop.hpp"

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
      f(i) = 1.0 / (1 + exp(beta_ * (x(i) - mu_)));
    ac(i) = acos(x(i));
  }
  for (unsigned i = 0; i < N; i++) {
    T sum = 0.0;
    for (unsigned j = 0; j < N; j++)
      sum += cos(i * ac(j)) * f(j);
    coef(i) = sum * 2.0 / N;
  }
  coef(0) *= 0.5;
  return coef;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_st_lcm()
{
  debug_message("Entered Simulation::calc_st_lcm\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_st_lcm = false;
    try {
      int dummy_variable;
      get_hdf5<
        int>(&dummy_variable, file, (char *)"/Calculation/STLCM/NumVectors");
      Global.calculate_st_lcm = true;
    } catch (H5::Exception &e) {
      debug_message("st_lcm: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_st_lcm = false;
#pragma omp critical
  local_calculate_st_lcm = Global.calculate_st_lcm;
#pragma omp barrier
  if (local_calculate_st_lcm) {
#pragma omp master
    std::cout << "Calculating Local Chern.\n";
#pragma omp barrier
    int number_samples;
    value_type energy;
    value_type beta;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<
        int>(&number_samples, file, (char *)"/Calculation/STLCM/NumVectors");
      get_hdf5<value_type>(&energy, file, (char *)"/Calculation/STLCM/Miu");
      get_hdf5<value_type>(&beta, file, (char *)"/Calculation/STLCM/Beta");
      file->close();
      delete file;
    }
    st_lcm(number_samples, beta, energy);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::st_lcm(
  const unsigned vectors_,
  const value_type beta_,
  const value_type energy_
)
{
  debug_message("Entered local chern\n");
  value_type energy_scale;
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
    file->close();
    delete file;
  }
#pragma omp barrier
  const value_type norm = 4 * M_PI / r.rLat.determinant();
  const value_type target = energy_ / energy_scale;
  const value_type beta = beta_ * energy_scale;
  const Eigen::Array<value_type, -1, 1> coefs = build_fermi(beta, target);

  KPM_Vector<T, D> phi(2, *this);
  KPM_Vector<T, D> ket(1, *this);
  KPM_Vector<T, D> bra(1, *this);
  std::vector<KPM_Vector<T, D> *> vectors;
  vectors.push_back(&phi);
  vectors.push_back(&ket);
  vectors.push_back(&bra);
  Eigen::Array<T, -1, -1> results(r.Sized, 2);
  results.setZero();
  Eigen::Array<T, -1, 1> prv(r.Sized);

  h.generate_disorder();
  for (unsigned vec = 0; vec < vectors_; ++vec) {
    for (auto &vec : vectors)
      vec->initiate_phases();
    phi.set_index(0);
    phi.initiate_vector();
    bra.v.col(0) = phi.v.col(0);
    phi.Exchange_Boundaries();
    ket.v.col(0).setZero();
    for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
      phi.cheb_iteration(n);
      ket.v.col(0) += phi.v.col(phi.get_index()) * coefs(n);
    }
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
    phi.set_index(0);
    phi.v.col(0) = ket.v.col(0);
    for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
      phi.cheb_iteration(n);
      ket.v.col(0) += phi.v.col(phi.get_index()) * coefs(n);
    }
    bra.empty_ghosts(0);
    const value_type weight = 1.0 / (vec + 1);
    const Eigen::Array<value_type, -1, 1> map =
      (bra.v.col(0) * ket.v.col(0)).imag();
    prv = results.col(0);
    results.col(0) += weight * (map - results.col(0));
    results.col(1) +=
      weight * ((map - prv) * (map - results.col(0)) - results.col(1));
  }
  results.col(1) = results.col(1).sqrt() / std::sqrt(vectors_);
  store_lcm(results);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_lcm(const Eigen::Array<T, -1, 1> &results_)
{
  debug_message("Entered store_ldos\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.lcm_map.resize(r.Sizet, 2);
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
      Global.lcm_map.row(global.index) = results_.row(local.index);
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/STLCM/Map");
    const std::string name(buffer);
    write_hdf5(Global.ldos_map, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_lcm\n");
}

#define INSTANTIATE_FERMI(type)                                                \
  template Eigen::Array<type, -1, 1> build_fermi(const type, const type);

INSTANTIATE_FERMI(float)
INSTANTIATE_FERMI(double)
INSTANTIATE_FERMI(long double)
#undef INSTANTIATE_FERMI

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
