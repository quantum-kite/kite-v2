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
void Simulation<T, D>::calc_custom_two_local()
{
  debug_message("Entered Simulation::calc_custom_two\n");
  std::string base_grp = "/Calculation/CustomTwoLocal/";
  std::string tmp = base_grp + "Vertex0/NumMoments";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    Global.calculate_custom_two_local = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, &file, tmp);
      Global.calculate_custom_two_local = true;
    } catch (H5::Exception &e) {
      debug_message("CustomTwoLocal: no need to calculate Custom Two Local.\n");
    }
    file.close();
  }
#pragma omp barrier
  bool local_calculate_custom_two_local = false;
#pragma omp critical
  local_calculate_custom_two_local = Global.calculate_custom_two_local;
#pragma omp barrier
  if (local_calculate_custom_two_local) {
#pragma omp master
    std::cout << "Calculating Custom Two Local\n";
#pragma omp barrier
    Eigen::Array<unsigned, -1, -1, Eigen::RowMajor> tt;
    Eigen::Array<unsigned, -1, -1> positions;
    std::array<unsigned, 2> number_moments;
    std::vector<Eigen::Array<T, -1, 1>> coefs(2);
    std::vector<std::vector<std::string>> stream(2);
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> orb_operators;
#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom Two: checking if we need to calculate Custom Two.\n"
      );
      std::string path;
      for (unsigned i = 0; i < 2; ++i) {
        const std::string base = base_grp + "Vertex" + std::to_string(i);
        path = base + "/NumCoefficients";
        int number_coefficients;
        get_hdf5<int>(&number_coefficients, &file, path);
        coefs.at(i).resize(number_coefficients);
        path = base + "/Coefficients";
        get_hdf5<T>(coefs[i].data(), &file, path);
        path = base + "/Operators";
        my_get_hdf5(stream[i], file, path);
        path = base + "/NumMoments";
        get_hdf5<unsigned>(&number_moments[i], &file, path);
        path = base_grp + "Positions";
        std::cout << path << std::endl;
        H5::DataSet dataset = file.openDataSet(path);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, nullptr);
        tt.resize(dims[0], dims[1]);
        get_hdf5<unsigned>(tt.data(), &file, path);
        positions = tt; // This transposes as tt is RowMajor
      }
      tmp = base_grp + "CustomOperators/";
      H5::Group grp;
      grp = file.openGroup(tmp);
      hsize_t n = grp.getNumObjs();
      for (hsize_t i = 0; i < n; ++i) {
        H5std_string memberName = grp.getObjnameByIdx(i);
        if (
          auto err =
            this->getMembers(grp, std::string(memberName), &orb_operators)
        ) {
        }
      }
      file.close();
    }
#pragma omp barrier
    custom_two_local(positions, number_moments, stream, coefs, orb_operators);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_two_local(
  const Eigen::Array<unsigned, -1, -1> &positions_,
  const std::array<unsigned, 2> &number_moments_,
  const std::vector<std::vector<std::string>> &stream_,
  const std::vector<Eigen::Array<T, -1, 1>> &coeffs_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_
)
{
  int size_gamma = 1;
  for (const auto &m : number_moments_) {
    if (m % 2 != 0) {
      std::cout << "The number of moments must be even\n";
      exit(1);
    }
    size_gamma *= m;
  }
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Eigen::Matrix<T, MEMORY, MEMORY> tmp;
  KPM_Vector<T, D> kpm_trc_0(1, *this);
  KPM_Vector<T, D> kpm_trc_1(2, *this);
  std::array<KPM_Vector<T, D> *, 2> vectors{&kpm_trc_0, &kpm_trc_1};
  KPM_Vector<T, D> kpm_trc_2(MEMORY, *this);
  KPM_Vector<T, D> kpm_trc_3(MEMORY, *this);

  h.generate_disorder();
  h.generate_twists();
  for (auto &vec : vectors)
    vec->initiate_phases();
  kpm_trc_2.initiate_phases();
  kpm_trc_3.initiate_phases();

  for (unsigned p = 0, P = positions_.rows(); p < P; ++p) {
    Eigen::Array<T, -1, -1> gamma(1, size_gamma);
    const Eigen::Array<unsigned, 1, D + 1> pos = positions_.row(p);
    if constexpr (D == 2)
      global.set({pos(1), pos(0), pos(2)});
    else
      global.set({pos(2), pos(1), pos(0), pos(3)});

    kpm_trc_0.build_site(global.index);
    kpm_trc_0.Exchange_Boundaries();
    kpm_trc_1.set_index(0);
    act_with_stream(stream_[0], operators_, coeffs_[0], vectors, 0);
    vectors[0] = &kpm_trc_1;
    vectors[1] = &kpm_trc_3;
    for (int n = 0, N = number_moments_[0]; n < N; n += MEMORY) {
      for (int i = n; i < n + MEMORY; ++i) {
        kpm_trc_1.cheb_iteration(i);
        act_with_stream(stream_[1], operators_, coeffs_[1], vectors, i % MEMORY);
        kpm_trc_3.empty_ghosts(i % MEMORY);
      }
      kpm_trc_2.set_index(0);
      kpm_trc_2.v.col(0) = kpm_trc_0.v.col(0);
      for (int m = 0, M = number_moments_[1]; m < M; m += MEMORY) {
        for (int i = m; i < m + MEMORY; i++)
          kpm_trc_2.cheb_iteration(i);
        tmp.setZero();
        for (std::size_t ii = 0; ii < r.Sized; ii += r.Ld[0])
          tmp += kpm_trc_2.v.block(ii, 0, r.Ld[0], MEMORY).adjoint() *
                 kpm_trc_3.v.block(ii, 0, r.Ld[0], MEMORY);
        for (unsigned j = 0; j < MEMORY; j++) {
          const unsigned idx = m + M * (n + j);
          for (unsigned i = 0; i < MEMORY; i++)
            gamma(idx + i) = tmp(i, j);
        }
      }
    }
    vectors[0] = &kpm_trc_0;
    vectors[1] = &kpm_trc_1;
    store_custom_two_local(p, number_moments_, pos, gamma);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_two_local(
  const unsigned p_,
  const std::array<unsigned, 2> &number_moments_,
  const Eigen::Array<unsigned, 1, D + 1> &pos_,
  const Eigen::Array<T, -1, -1> &gamma_
)
{
  const std::string base_grp =
    "/Calculation/CustomTwoLocal/p_" + std::to_string(p_) + "/";
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDWR);
    H5::Group grp;
    try {
      grp = file.openGroup(base_grp);
    } catch (...) {
      grp = file.createGroup(base_grp);
    }
  }
#pragma omp barrier
  Eigen::Array<T, -1, -1> general_gamma = Eigen::Map<const Eigen::Array<
    T, -1, -1>>(gamma_.data(), number_moments_[1], number_moments_[0]);
#pragma omp master
  {
    Global.general_gamma.resize(general_gamma.rows(), general_gamma.cols());
    Global.general_gamma.setZero();
  }
#pragma omp barrier
#pragma omp critical
  {
    Global.general_gamma.matrix() += general_gamma.matrix();
  }
#pragma omp barrier
#pragma omp master
  {
    const std::string np = base_grp + "Pos";
    H5::H5File file(name, H5F_ACC_RDWR);
    write_hdf5(Eigen::Array<unsigned, -1, -1>(pos_), &file, np);
    const std::string ng = base_grp + "Gamma";
    write_hdf5(Global.general_gamma, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_gamma\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
