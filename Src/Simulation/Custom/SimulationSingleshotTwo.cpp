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
T gauss_first(const T n_, const T mu_, const T sigma_)
{
  const T numerator = n_ * mu_ * sigma_ * sigma_ *
                      (n_ * n_ * sigma_ * sigma_ / (1 - mu_ * mu_) - 3);
  const T denominator = std::pow(1 - mu_ * mu_, -1.5);
  return numerator * denominator;
}

template <typename T>
T gauss_second(const T n_, const T mu_, const T sigma_)
{
  const T term_1 = 7 * mu_ * mu_ - 4;
  const T tmp = 1 - mu_ * mu_;
  const T term_2 = 3 - 6 * n_ * n_ * sigma_ * sigma_ / tmp +
                   std::pow(n_ * sigma_, 4) / (tmp * tmp);
  const T denominator = 1 / (24 * tmp);
  return sigma_ * sigma_ * term_1 * term_2 * denominator;
}

template <typename T>
Eigen::Array<T, -1, 1> build_gaussian(
  const T energy_,
  const T sigma_,
  const T escale_,
  const unsigned Ncheb_
)
{
  const T sigma = sigma_ / escale_;
  const T energy = energy_ / escale_;
  Eigen::Array<T, -1, 1> coefs(Ncheb_);
  coefs(0) = 1 - gauss_second(static_cast<T>(0), energy, sigma);
  for (unsigned n = 1; n < Ncheb_; ++n) {
    const T gaussian =
      std::exp(-0.5 * n * n * sigma * sigma / (1 - energy * energy));
    const T cossine = std::cos(n * std::acos(energy));
    const T sine = std::sin(n * std::acos(energy));
    coefs(n) = 2 * gaussian *
               (cossine * (1 - gauss_second(static_cast<T>(n), energy, sigma)) -
                0.5 * sine * gauss_first(static_cast<T>(n), energy, sigma));
  }
  const T prefactor = 1 / (M_PI * std::sqrt(1 - energy * energy));
  coefs *= prefactor;
  return coefs / escale_;
}

template <typename T>
Eigen::Array<std::complex<T>, -1, 1>
build_dgreen(const std::complex<T> z_, const T escale_, const unsigned Ncheb_)
{
  using cplx = std::complex<T>;
  constexpr cplx I{0.0, 1.0};
  Eigen::Array<cplx, -1, 1> coef(Ncheb_);

  const cplx z = z_ / escale_;
  const cplx sq = static_cast<T>(1.0) - z * z;
  const cplx sqr = std::sqrt(sq);

  const cplx diff = z - I * sqr;

  cplx current_power = 1.0;
  coef(0) = -I * z * current_power / (sq * sqr);
  current_power *= diff;
  for (unsigned n = 1; n < Ncheb_; ++n) {
    coef(n) = static_cast<T>(2.0) * (static_cast<T>(n) - I * z / sqr) *
              current_power / sq;
    current_power *= diff;
  }
  coef /= (escale_ * escale_);
  return coef;
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
    value_type energy_scale;
    std::array<Eigen::Array<value_type, -1, 1>, 3> targets;
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

      std::array<std::string, 3> fields{"Gamma", "Sigma", "Energies"};
      for (unsigned i = 0; i < fields.size(); ++i) {
        tmp = base_grp + fields[i];
        H5::DataSet dataset = file->openDataSet(tmp);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1] = {0};
        dataspace.getSimpleExtentDims(dims, nullptr);
        targets[i].resize(dims[0]);
        get_hdf5<value_type>(targets[i].data(), file, tmp);
      }
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
      number_samples, number_vectors, targets, number_moments, stream, coefs,
      orb_operators
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_ss_two(
  const int samples_,
  const int vectors_,
  const std::array<Eigen::Array<value_type, -1, 1>, 3> &targets_,
  const std::vector<int> &num_pol_,
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
    const std::array<Eigen::Index, 3>
      dims{targets_[0].size(), targets_[1].size(), targets_[2].size()};
    std::array<Eigen::Array<T, -1, -1>, 2> coefs;
    for (unsigned i = 0; i < 2; ++i)
      coefs[i].resize(num_pol_[i], dims[i] * dims[2]);

    for (unsigned i = 0; i < dims[2]; ++i) {
      unsigned count = 0;
      const value_type en = targets_[2](i);
      for (const auto &gamma : targets_[0]) {
        const unsigned idx = count + i * dims[0];
        const T z = en + I * gamma;
        coefs[0].col(idx) = build_dgreen(z, energy_scale, num_pol_[0]);
        ++count;
      }
      count = 0;
      for (const auto &sigma : targets_[1]) {
        const unsigned idx = count + i * dims[1];
        coefs[1].col(idx) =
          build_gaussian(en, sigma, energy_scale, num_pol_[1]);
        ++count;
      }
    }
    std::array<Eigen::Array<T, -1, 1>, 2> ket =
      {Eigen::Array<T, -1, 1>(r.Sized), Eigen::Array<T, -1, 1>(r.Sized)};
    Eigen::Array<T, -1, 1> bra(r.Sized);
    Eigen::Array<T, -1, 1> psi(r.Sized);
    Eigen::Array<T, -1, -1> result(dims[2], dims[0] * dims[1]);
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

        for (unsigned j = 0, J = dims[0] * dims[1]; j < J; ++j) {
          const unsigned idx_gamma = j % dims[0];
          const unsigned idx_sigma = j / dims[0];
          for (unsigned i = 0; i < dims[2]; ++i) {
            const std::array<Eigen::Index, 2>
              broad_id{idx_gamma + i * dims[0], idx_sigma + i * dims[1]};
            for (unsigned p = 0; p < 2; ++p) {
              vectors[3].set_index(0);
              vectors[3].v.col(0) = psi;
              ket[p].setZero();
              const unsigned id_0 = p % 2;
              const unsigned id_1 = (p + 1) % 2;
              for (unsigned n = 0, N = num_pol_[id_0]; n < N; ++n) {
                vectors[3].cheb_iteration(n);
                const unsigned idx_read = vectors[3].get_index();
                const T coef = coefs[id_0](n, broad_id[id_0]);
                ket[p] += coef * vectors[3].v.col(idx_read).array();
              }
              vectors[0].v.col(0) = ket[p];
              act_with_stream(
                stream_[1], operators_, coeffs_[1], vector_ptrs, 0
              );
              vectors[3].set_index(0);
              vectors[3].Exchange_Boundaries();
              ket[p].setZero();
              for (unsigned n = 0, N = num_pol_[id_1]; n < N; ++n) {
                vectors[3].cheb_iteration(n);
                const unsigned idx_read = vectors[3].get_index();
                const T coef = std::conj(coefs[id_1](n, broad_id[id_1]));
                ket[p] += coef * vectors[3].v.col(idx_read).array();
              }
            }
            T res = (bra * ket[0] + (bra * ket[1]).conjugate()).sum();
            res *= 0.5;
            result(i, j) += (res - result(i, j)) * weight;
          }
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

#define INSTANTIATE_GAUSS(type)                                                \
  template type gauss_first<type>(const type, const type, const type);         \
  template type gauss_second<type>(const type, const type, const type);        \
  template Eigen::Array<type, -1, 1>                                           \
  build_gaussian(const type, const type, const type, const unsigned);          \
  template Eigen::Array<std::complex<type>, -1, 1>                             \
  build_dgreen(const std::complex<type>, const type, const unsigned);

INSTANTIATE_GAUSS(float)
INSTANTIATE_GAUSS(double)
INSTANTIATE_GAUSS(long double)
#undef INSTANTIATE_GAUSS
