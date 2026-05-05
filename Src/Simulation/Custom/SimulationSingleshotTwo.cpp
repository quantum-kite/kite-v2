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

template <typename T>
std::complex<T> green(const unsigned n_, const std::complex<T> z_)
{
  using cplx = std::complex<T>;
  constexpr cplx I{0.0, 1.0};
  const T norm = -2.0 / (1.0 + static_cast<unsigned>(n_ == 0));
  const cplx sqr = std::sqrt(static_cast<T>(1.0) - z_ * z_);
  const cplx diff = z_ - I * sqr;
  const cplx power = std::pow(diff, static_cast<int>(n_));
  return I * norm * power / sqr;
}

template <typename T>
std::complex<T> dgreen(const unsigned n_, const std::complex<T> z_)
{
  using cplx = std::complex<T>;
  constexpr cplx I{0.0, 1.0};
  const T norm = 2.0 / (1.0 + static_cast<unsigned>(n_ == 0));
  const cplx sq = static_cast<T>(1.0) - z_ * z_;
  const cplx sqr = std::sqrt(sq);
  const cplx diff = z_ - I * sqr;
  const cplx power = std::pow(diff, static_cast<int>(n_));
  return norm * (static_cast<T>(n_) - I * z_ / sqr) * power / sq;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_custom_ss_two()
{
  debug_message("Entered Simulation::calc_custom_single_two\n");
  std::string base_grp = "/Calculation/CustomSingleTwo/";
  std::string tmp = base_grp + "NumMoments";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_custom_ss_two = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, file, tmp);
      Global.calculate_custom_ss_two = true;
    } catch (H5::Exception &e) {
      debug_message(
        "CustomSingleTwo: no need to calculate Custom Singleshot Two.\n"
      );
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_custom_ss_two = false;
#pragma omp critical
  local_calculate_custom_ss_two = Global.calculate_custom_ss_two;
#pragma omp barrier
  if (local_calculate_custom_ss_two) {
#pragma omp master
    std::cout << "Calculating Custom Singleshot Two\n";
#pragma omp barrier
    int number_samples;
    int number_vectors;
    Eigen::Array<value_type, -1, 1> energies;
    Eigen::Array<value_type, -1, 1> gammas;
    std::vector<int> number_moments(2);
    std::vector<Eigen::Array<T, -1, 1>> coefs(2);
    std::vector<std::vector<std::string>> stream(2);
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> orb_operators;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom Single Two: checking if we need to calculate Custom Single "
        "Two.\n"
      );
      tmp = base_grp + "NumDisorder";
      get_hdf5<int>(&number_samples, file, tmp);
      tmp = base_grp + "NumVectors";
      get_hdf5<int>(&number_vectors, file, tmp);

      tmp = base_grp + "Energies";
      H5::DataSet dataset = file->openDataSet(tmp);
      H5::DataSpace dataspace = dataset.getSpace();
      hsize_t dims[1] = {0};
      dataspace.getSimpleExtentDims(dims, nullptr);
      energies.resize(dims[0]);
      get_hdf5<value_type>(energies.data(), file, tmp);

      tmp = base_grp + "Gamma";
      dataset = file->openDataSet(tmp);
      dataspace = dataset.getSpace();
      dataspace.getSimpleExtentDims(dims, nullptr);
      gammas.resize(dims[0]);
      get_hdf5<value_type>(gammas.data(), file, tmp);

      std::string path = base_grp + "/NumMoments";
      for (unsigned i = 0; i < 2; ++i)
        get_hdf5<int>(&number_moments[i], file, path);

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
    custom_ss_two(
      number_samples, number_vectors, energies, gammas, number_moments, stream,
      coefs, orb_operators
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_ss_two(
  const int samples_,
  const int vectors_,
  const Eigen::Array<value_type, -1, 1> &energies_,
  const Eigen::Array<value_type, -1, 1> &gammas_,
  const std::vector<int> &number_moments_,
  const std::vector<std::vector<std::string>> &stream_,
  const std::vector<Eigen::Array<T, -1, 1>> &coeffs_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_
)
{
  if constexpr (is_tt<std::complex, T>::value) {
    constexpr T I{0.0, 1.0};
    value_type energy_scale;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
      file->close();
      delete file;
    }
#pragma omp barrier
    const unsigned en_size = energies_.size();
    const unsigned col_num = en_size * gammas_.size();
    std::vector<std::array<Eigen::Array<T, -1, 1>, 2>> coefs(col_num);
    for (unsigned i = 0; i < col_num; ++i) {
      const unsigned idx_e = i % en_size;
      const unsigned idx_g = i / en_size;
      const value_type energy = energies_(idx_e) / energy_scale;
      const value_type gamma = gammas_(idx_g) / energy_scale;
      const T zz = energy + I * gamma;
      coefs[i][0] = Eigen::Array<T, -1, 1>::Zero(number_moments_[0]);
      coefs[i][1] = Eigen::Array<T, -1, 1>::Zero(number_moments_[1]);

      for (unsigned a = 0, A = number_moments_[0]; a < A; ++a)
        coefs[i][0](a) = -green<value_type>(a, zz).imag() / M_PI;
      for (unsigned a = 0, A = number_moments_[1]; a < A; ++a)
        coefs[i][1](a) = dgreen<value_type>(a, zz);
    }
    std::array<Eigen::Array<T, -1, 1>, 2> ket =
      {Eigen::Array<T, -1, 1>(r.Sized), Eigen::Array<T, -1, 1>(r.Sized)};
    Eigen::Array<T, -1, 1> bra(r.Sized);
    Eigen::Array<T, -1, 1> psi(r.Sized);
    Eigen::Array<T, -1, -1> result(en_size, gammas_.size());
    result.setZero();
    std::array<KPM_Vector<T, D>, 4> vectors =
      {KPM_Vector<T, D>(1, *this), KPM_Vector<T, D>(1, *this),
       KPM_Vector<T, D>(1, *this), KPM_Vector<T, D>(2, *this)};
    std::vector<KPM_Vector<T, D> *> vector_ptrs =
      {&vectors[0], &vectors[1], &vectors[2], &vectors[3]};

    unsigned average = 0;
    for (int disorder = 0; disorder < samples_; ++disorder) {
      h.generate_disorder();
      for (int vec = 0; vec < vectors_; ++vec) {
        h.generate_twists();
        vectors[0].initiate_vector();
        bra = vectors[0].v.col(0).conjugate();
        vectors[0].Exchange_Boundaries();
        for (auto &vec : vectors)
          vec.initiate_phases();
        vectors[3].set_index(0);

        // Act with stream[0]
        act_with_stream(stream_[0], operators_, coeffs_[0], vector_ptrs, 0);
        psi = vectors[3].v.col(0);

        const value_type weight = 1.0 / static_cast<value_type>(average + 1);
        for (unsigned i = 0; i < col_num; ++i) {
          const unsigned idx_e = i % en_size;
          const unsigned idx_g = i / en_size;
          const std::array<Eigen::Array<T, -1, 1>, 2> &cfs = coefs[i];

          for (unsigned j = 0; j < 2; ++j) {
            vectors[3].set_index(0);
            vectors[3].v.col(0) = psi;
            ket[j].setZero();
            const unsigned c0_idx = (j + 1) % 2;
            const unsigned c1_idx = j % 2;
            for (unsigned n = 0, N = number_moments_[c0_idx]; n < N; ++n) {
              vectors[3].cheb_iteration(n);
              const T coef = cfs[c0_idx](n);
              ket[j] += coef * vectors[3].v.col(vectors[3].get_index()).array();
            }
            vectors[0].v.col(0) = ket[j];
            act_with_stream(stream_[1], operators_, coeffs_[1], vector_ptrs, 0);

            vectors[3].set_index(0);
            vectors[3].Exchange_Boundaries();
            ket[j].setZero();
            for (unsigned n = 0, N = number_moments_[c1_idx]; n < N; ++n) {
              vectors[3].cheb_iteration(n);
              const T coef = std::conj(cfs[c1_idx](n));
              ket[j] += coef * vectors[3].v.col(vectors[3].get_index()).array();
            }
          }
          const T res = (bra * ket[0]).sum() + (bra * ket[1]).conjugate().sum();
          result(idx_e, idx_g) += weight * (static_cast<value_type>(0.5) * res -
                                            result(idx_e, idx_g));
        }
        ++average;
      }
    }
    store_custom_ss_two(result, average);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_ss_two(
  const Eigen::Array<T, -1, -1> &result_,
  const unsigned avr_
)
{
  debug_message("Entered store_gamma\n");
#pragma omp master
  {
    Global.singleshot_cond.resize(result_.rows(), result_.cols());
    Global.singleshot_cond.setZero();
  }
#pragma omp barrier
#pragma omp critical
  Global.singleshot_cond += result_;
#pragma omp barrier
#pragma omp master
  {
    const std::array<std::string, 2> names{
      "/Calculation/CustomSingleTwo/Result", "/Calculation/CustomSingleTwo/Avr"
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
        write_hdf5(Global.singleshot_cond, &file, n);
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

#define INSTANTIATE_GREEN(type)                                                \
  template std::complex<type>                                                  \
  dgreen(const unsigned, const std::complex<type>);                            \
  template std::complex<type> green(const unsigned, const std::complex<type>);

INSTANTIATE_GREEN(float)
INSTANTIATE_GREEN(double)
INSTANTIATE_GREEN(long double)
#undef INSTANTIATE_GREEN
