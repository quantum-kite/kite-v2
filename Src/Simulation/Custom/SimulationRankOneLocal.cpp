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

template <typename T, unsigned D>
void Simulation<T, D>::calc_custom_one_local()
{
  debug_message("Entered Simulation::calc_custom_one_local\n");
  std::string base_grp = "/Calculation/CustomOneLocal/";
  std::string tmp = base_grp + "NumMoments";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_custom_density = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, file, tmp);
      Global.calculate_custom_density = true;
    } catch (H5::Exception &e) {
      debug_message(
        "Custom One - Local: no need to calculate Custom One - Local.\n"
      );
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_custom_density = false;
#pragma omp critical
  local_calculate_custom_density = Global.calculate_custom_density;
#pragma omp barrier
  if (local_calculate_custom_density) {
#pragma omp master
    std::cout << "Calculating Custom One - Local\n";
#pragma omp barrier
    int number_samples;
    int number_moments;
    int number_local_sites;
    std::vector<int> orbitals;
    std::vector<int> positions;
    Eigen::Array<T, -1, 1> coefs;
    std::vector<std::string> stream;
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>>
      operator_collection;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom One - Local: checking if we need to calculate Custom One - "
        "Local.\n"
      );
      tmp = base_grp + "NumDisorder";
      get_hdf5<int>(&number_samples, file, tmp);
      tmp = base_grp + "NumMoments";
      get_hdf5<int>(&number_moments, file, tmp);

      tmp = base_grp + "Orbitals";
      H5::DataSet orb_dset = file->openDataSet(tmp);
      H5::DataSpace orb_dspace = orb_dset.getSpace();
      hsize_t orb_dims[2];
      orb_dspace.getSimpleExtentDims(orb_dims, nullptr);
      orbitals.resize(orb_dims[0]);
      get_hdf5(orbitals.data(), file, tmp);

      tmp = base_grp + "FixPosition";
      H5::DataSet pos_dset = file->openDataSet(tmp);
      H5::DataSpace pos_dspace = pos_dset.getSpace();
      hsize_t pos_dims[2];
      pos_dspace.getSimpleExtentDims(pos_dims, nullptr);
      positions.resize(pos_dims[0]);
      get_hdf5(positions.data(), file, tmp);
      number_local_sites = orbitals.size();

      const std::string base = base_grp + "Vertex0";
      std::string path = base + "/NumCoefficients";
      int number_coefficients;
      get_hdf5<int>(&number_coefficients, file, path);
      coefs.resize(number_coefficients);
      path = base + "/Coefficients";
      get_hdf5<T>(coefs.data(), file, path);
      path = base + "/Operators";
      my_get_hdf5(stream, my_file, path);

      tmp = base_grp + "CustomOperators/";
      H5::Group grp;
      grp = file->openGroup(tmp);
      hsize_t n = grp.getNumObjs();
      for (hsize_t i = 0; i < n; ++i) {
        H5std_string memberName = grp.getObjnameByIdx(i);
        if (
          auto err =
            this->getMembers(grp, std::string(memberName), &operator_collection)
        ) {
        }
      }
      file->close();
      delete file;
    }
#pragma omp barrier
    custom_one_local(
      number_moments, number_samples, number_local_sites, orbitals, positions,
      coefs, operator_collection, stream
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_one_local(
  const int number_moments_,
  const int number_samples_,
  const int number_local_sites_,
  const std::vector<int> &orbitals_,
  const std::vector<int> &positions_,
  const Eigen::Array<T, -1, 1> &coeffs_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_,
  const std::vector<std::string> &stream_
)
{
  Eigen::Matrix<T, -1, -1> gamma(number_moments_, number_local_sites_);
  gamma.setZero();
  Eigen::Matrix<T, 1, 2> tmp;
  tmp.setZero();
  Eigen::Matrix<T, -1, 1> gamma_site(number_moments_);
  KPM_Vector<T, D> kpm_trc_0(1, *this);
  KPM_Vector<T, D> kpm_trc_1(2, *this);
  std::array<KPM_Vector<T, D> *, 2> vectors{&kpm_trc_0, &kpm_trc_1};
  h.generate_twists();
  for (int disorder = 0; disorder < number_samples_; ++disorder) {
    h.generate_disorder();
    for (int site = 0; site < number_local_sites_; ++site) {
      gamma_site.setZero();
      const std::size_t global_index =
        positions_[site] * r.Ld[0] + orbitals_[site];
      kpm_trc_0.build_site(global_index);
      kpm_trc_0.Exchange_Boundaries();
      for (auto &vec : vectors)
        vec->initiate_phases();
      kpm_trc_1.set_index(0);
      act_with_stream(stream_, operators_, coeffs_, vectors, 0);
      kpm_trc_0.empty_ghosts(0);
      for (int m = 0; m < number_moments_; m += 2) {
        kpm_trc_1.cheb_iteration(m);
        kpm_trc_1.cheb_iteration(m + 1);
        tmp.setZero();
        for (std::size_t ii = 0; ii < r.Sized; ii += r.Ld[0]) {
          tmp += kpm_trc_0.v.block(ii, 0, r.Ld[0], 1).adjoint() *
                 kpm_trc_1.v.block(ii, 0, r.Ld[0], 2);
        }
        gamma_site.segment(m, 2) = tmp;
      }
      gamma.col(site) += (gamma_site - gamma.col(site)) / (disorder + 1);
    }
  }
  store_custom_one_local(gamma, 0);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_one_local(
  const Eigen::Matrix<T, -1, -1> &gamma_,
  const unsigned avr_
)
{
#pragma omp master
  Global.general_gamma =
    Eigen::Array<T, -1, -1>::Zero(gamma_.rows(), gamma_.cols());
#pragma omp barrier
#pragma omp critical
  Global.general_gamma.matrix() += 0.5 * (gamma_ + gamma_.conjugate());
#pragma omp barrier
#pragma omp master
  {
    const std::array<std::string, 2> names{
      "/Calculation/CustomOneLocal/Gamma", "/Calculation/CustomOneLocal/Avr"
    };
    Eigen::Array<value_type, -1, -1> avr_arr;
    avr_arr.resize(1, 1);
    avr_arr(0) = avr_;
    unsigned count = 0;
    for (const auto &n : names) {
      H5::H5File file(name, H5F_ACC_RDWR);
      if (H5Lexists(file.getId(), n.c_str(), H5P_DEFAULT))
        file.unlink(n);
      if (count == 0)
        write_hdf5(Global.general_gamma, &file, n);
      else
        write_hdf5(avr_arr, &file, n);
      ++count;
    }
  }
#pragma omp barrier
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
