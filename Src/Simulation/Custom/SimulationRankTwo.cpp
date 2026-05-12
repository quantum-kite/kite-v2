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
void Simulation<T, D>::calc_custom_two()
{
  debug_message("Entered Simulation::calc_custom_two\n");
  std::string base_grp = "/Calculation/CustomTwo/";
  std::string tmp = base_grp + "NumDisorder";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_custom_conddc = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, file, tmp);
      Global.calculate_custom_conddc = true;
    } catch (H5::Exception &e) {
      debug_message("CustomTwo: no need to calculate Custom Two.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_custom_conddc = false;
#pragma omp critical
  local_calculate_custom_conddc = Global.calculate_custom_conddc;
#pragma omp barrier
  if (local_calculate_custom_conddc) {
#pragma omp master
    std::cout << "Calculating Custom Two\n";
#pragma omp barrier
    int number_samples;
    int number_vectors;
    std::vector<int> number_moments(2);
    std::vector<Eigen::Array<T, -1, 1>> coefs(2);
    std::vector<std::vector<std::string>> stream(2);
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> orb_operators;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom Two: checking if we need to calculate Custom Two.\n"
      );
      tmp = base_grp + "NumDisorder";
      get_hdf5<int>(&number_samples, file, tmp);
      tmp = base_grp + "NumVectors";
      get_hdf5<int>(&number_vectors, file, tmp);

      std::string path;
      for (unsigned i = 0; i < 2; ++i) {
        const std::string base = base_grp + "Vertex" + std::to_string(i);
        path = base + "/NumCoefficients";
        int number_coefficients;
        get_hdf5<int>(&number_coefficients, file, path);
        coefs.at(i).resize(number_coefficients);
        path = base + "/Coefficients";
        get_hdf5<T>(coefs.at(i).data(), file, path);
        path = base + "/Operators";
        my_get_hdf5(stream.at(i), my_file, path);
        path = base + "/NumMoments";
        get_hdf5<int>(&number_moments[i], file, path);
      }
      tmp = base_grp + "CustomOperators/";
      H5::Group grp;
      grp = file->openGroup(tmp);
      hsize_t n = grp.getNumObjs();
      for (hsize_t i = 0; i < n; ++i) {
        H5std_string memberName = grp.getObjnameByIdx(i);
        if (
          auto err =
            this->getMembers(grp, std::string(memberName), &orb_operators)
        ) {
        }
      }
      file->close();
      delete file;
    }
#pragma omp barrier
    custom_two(
      number_samples, number_vectors, number_moments, stream, coefs,
      orb_operators
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_two(
  const int samples_,
  const int random_vectors_,
  const std::vector<int> &number_moments_,
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
  Eigen::Matrix<T, MEMORY, MEMORY> tmp;
  Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1>::Zero(1, size_gamma);

  KPM_Vector<T, D> kpm_trc_0(1, *this);
  KPM_Vector<T, D> kpm_trc_1(2, *this);
  std::array<KPM_Vector<T, D> *, 2> vectors{&kpm_trc_0, &kpm_trc_1};
  KPM_Vector<T, D> kpm_trc_2(MEMORY, *this);
  KPM_Vector<T, D> kpm_trc_3(MEMORY, *this);
  unsigned average = 0;

  for (int rand_v = 0; rand_v < random_vectors_; ++rand_v) {
    h.generate_twists();
    kpm_trc_0.initiate_vector();
    kpm_trc_0.Exchange_Boundaries();
    for (int disorder = 0; disorder < samples_; ++disorder) {
      h.generate_disorder();
      for (unsigned i = 0, I = vectors.size(); i < I; ++i)
        vectors[i]->initiate_phases();
      kpm_trc_2.initiate_phases();
      kpm_trc_3.initiate_phases();
      kpm_trc_1.set_index(0);

      act_with_stream(stream_[0], operators_, coeffs_[0], vectors, 0);
      vectors[0] = &kpm_trc_1;
      vectors[1] = &kpm_trc_3;
      for (int n = 0, N = number_moments_[0]; n < N; n += MEMORY) {
        for (int i = n; i < n + MEMORY; ++i) {
          kpm_trc_1.cheb_iteration(i);
          act_with_stream(
            stream_[1], operators_, coeffs_[1], vectors, i % MEMORY
          );
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
          for (int j = 0; j < MEMORY; j++)
            for (int i = 0; i < MEMORY; i++) {
              const int ind = m + i + M * (n + j);
              gamma(ind) += (tmp(i, j) - gamma(ind)) / value_type(average + 1);
            }
        }
      }
      ++average;
      vectors[0] = &kpm_trc_0;
      vectors[1] = &kpm_trc_1;
    }
  }
  store_custom_two(gamma, average, number_moments_);
}

// Method for exact trace calculation
// template <typename T, unsigned D>
// void Simulation<T, D>::custom_two(
//   const int samples_,
//   const int random_vectors_,
//   const std::vector<int> &number_moments_,
//   const std::vector<std::vector<std::string>> &stream_,
//   const std::vector<Eigen::Array<T, -1, 1>> &coeffs_,
//   const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_
// )
// {
//   int size_gamma = 1;
//   for (const auto &m : number_moments_) {
//     if (m % 2 != 0) {
//       std::cout << "The number of moments must be even\n";
//       exit(1);
//     }
//     size_gamma *= m;
//   }
//   Eigen::Matrix<T, MEMORY, MEMORY> tmp;
//   Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1>::Zero(1,
//   size_gamma);

//   KPM_Vector<T, D> kpm_trc_0(1, *this);
//   KPM_Vector<T, D> kpm_trc_1(2, *this);
//   KPM_Vector<T, D> kpm_aux_0(1, *this);
//   KPM_Vector<T, D> kpm_aux_1(1, *this);
//   std::vector<KPM_Vector<T, D> *> vectors;
//   vectors.push_back(&kpm_trc_0);
//   vectors.push_back(&kpm_aux_0);
//   vectors.push_back(&kpm_aux_1);
//   vectors.push_back(&kpm_trc_1);
//   for (unsigned i = 0, I = vectors.size(); i < I; ++i)
//     vectors.at(i)->initiate_phases();
//   KPM_Vector<T, D> kpm_trc_2(MEMORY, *this);
//   KPM_Vector<T, D> kpm_trc_3(MEMORY, *this);
//   kpm_trc_2.initiate_phases();
//   kpm_trc_3.initiate_phases();
//   h.generate_twists();
//   h.generate_disorder();
//   for (std::size_t pos = 0, Pos = r.Sizet; pos < Pos; ++pos) {
//     kpm_trc_0.build_site(pos);
//     kpm_trc_0.Exchange_Boundaries();
//     kpm_trc_1.set_index(0);

//     act_with_stream(stream_.at(0), operators_, coeffs_.at(0), vectors, 0);
//     vectors.at(0) = &kpm_trc_1;
//     vectors.at(3) = &kpm_trc_3;
//     for (int n = 0, N = number_moments_.at(0); n < N; n += MEMORY) {
//       for (int i = n; i < n + MEMORY; ++i) {
//         kpm_trc_1.cheb_iteration(i);
//         act_with_stream(
//           stream_.at(1), operators_, coeffs_.at(1), vectors, i % MEMORY
//         );
//         kpm_trc_3.empty_ghosts(i % MEMORY);
//       }
//       kpm_trc_2.set_index(0);
//       kpm_trc_2.v.col(0) = kpm_trc_0.v.col(0);
//       for (int m = 0, M = number_moments_.at(1); m < M; m += MEMORY) {
//         for (int i = m; i < m + MEMORY; i++)
//           kpm_trc_2.cheb_iteration(i);
//         tmp.setZero();
//         for (std::size_t ii = 0; ii < r.Sized; ii += r.Ld[0])
//           tmp += kpm_trc_2.v.block(ii, 0, r.Ld[0], MEMORY).adjoint() *
//                  kpm_trc_3.v.block(ii, 0, r.Ld[0], MEMORY);
//         for (int j = 0; j < MEMORY; j++)
//           for (int i = 0; i < MEMORY; i++)
//             gamma(m + i + M * (n + j)) += tmp(i, j);
//       }
//     }
//     vectors.at(0) = &kpm_trc_0;
//     vectors.at(3) = &kpm_trc_1;
//   }
//   store_custom_two(gamma, 1, number_moments_);
// }

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_two(
  Eigen::Array<T, -1, -1> &gamma_,
  const unsigned avr_,
  const std::vector<int> &number_moments_
)
{
  debug_message("Entered store_gamma\n");
  Eigen::Array<T, -1, -1> general_gamma = Eigen::Map<Eigen::Array<
    T, -1, -1>>(gamma_.data(), number_moments_[1], number_moments_[0]);
#pragma omp master
  {
    Global.general_gamma.resize(general_gamma.rows(), general_gamma.cols());
    Global.general_gamma.setZero();
  }
#pragma omp barrier
#pragma omp critical
  {
    if (number_moments_[0] == number_moments_[1])
      Global.general_gamma.matrix() +=
        0.5 * (general_gamma.matrix() + general_gamma.matrix().adjoint());
    else
      Global.general_gamma.matrix() += general_gamma.matrix();
  }
#pragma omp barrier
#pragma omp master
  {
    const std::array<std::string, 2>
      names{"/Calculation/CustomTwo/Gamma", "/Calculation/CustomTwo/Avr"};
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
  debug_message("Left store_gamma\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
