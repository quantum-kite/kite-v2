/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Operator-weighted generalization of the stochastic/Markov LDOS map method
// (SimulationLDoS.cpp): instead of the plain per-site density
// map_r = factor * |<r|ket>|^2, computes, for each user-registered on-site
// operator O_r, map_r^O = factor * Re(<ket|_r O_r |ket>_r), using only the
// propagated ket (never the raw random seed) at both orbital indices -- the
// unbiased estimator derived for this feature. Reduces exactly to the plain
// formula when O_r is a single-orbital projector. Kept in its own
// translation unit (mirrors custom_one/custom_two) so the plain
// calc_ldos/ldos/store_ldos path is untouched.

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
#include "Coefficients.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_ldos_operators()
{
  debug_message("Entered Simulation::calc_ldos_operators\n");
#pragma omp barrier

  std::vector<std::string> operator_labels;
  std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> operator_collection;
  int vectors = 0, coef_id = 0;
  value_type energy = 0, sigma = 0;
  bool has_operators = false;

#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    try {
      std::string tmp = "/Calculation/ldos_map/Operators";
      my_get_hdf5(operator_labels, *file, tmp);

      get_hdf5<int>(&vectors, file, (char *) "/Calculation/ldos_map/NumVectors");
      get_hdf5<value_type>(&energy, file, (char *) "/Calculation/ldos_map/Energy");
      get_hdf5<value_type>(&sigma, file, (char *) "/Calculation/ldos_map/Sigma");
      get_hdf5<int>(&coef_id, file, (char *) "/Calculation/ldos_map/Coef_ID");

      H5::Group grp = file->openGroup("/Calculation/ldos_map/CustomOperators/");
      for (const auto &label : operator_labels) {
        if (auto err = this->getMembers(grp, label, &operator_collection)) {
          std::cerr << "getMembers failed for " << label << "\n";
        }
      }
      has_operators = !operator_labels.empty();
    } catch (H5::Exception &e) {
      debug_message("ldos_map: no operator-weighted map requested.\n");
    }
    file->close();
    delete file;
    Global.calculate_ldos_map_operators = has_operators;
  }
#pragma omp barrier

  bool local_has_operators;
#pragma omp critical
  local_has_operators = Global.calculate_ldos_map_operators;
#pragma omp barrier

  if (local_has_operators) {
#pragma omp master
    std::cout << "Calculating operator-weighted LDoS map.\n";
#pragma omp barrier
    ldos_operators(vectors, energy, sigma, coef_id, operator_labels, operator_collection);
  }
  debug_message("Left Simulation::calc_ldos_operators\n");
}

template <typename T, unsigned D>
void Simulation<T, D>::ldos_operators(
  const int vectors_,
  const value_type energy_,
  const value_type sigma_,
  const int coef_id_,
  const std::vector<std::string> &operator_labels,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operator_collection
)
{
  debug_message("Entered ldos_operators\n");
  if constexpr (is_tt<std::complex, T>::value) {
    const unsigned NOperators = static_cast<unsigned>(operator_labels.size());

    value_type energy_scale;
    value_type energy_shift;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(&energy_scale, file, (char *) "/EnergyScale");
      get_hdf5<value_type>(&energy_shift, file, (char *) "/EnergyShift");
      file->close();
      delete file;
    }
#pragma omp barrier

    const value_type target = (energy_ - energy_shift) / energy_scale;
    const value_type sigma = sigma_ / energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type factor =
      (coef_id_) ? 1.0 : std::sqrt(8 * M_PI) * sigma / energy_scale;
    const value_type fwhm = (coef_id_) ? sigma : std::sqrt(2) * sigma;
    const Eigen::Array<value_type, -1, 1> coefs =
      (coef_id_) ? Coefficients::build_window<value_type>(target, fwhm)
                 : Coefficients::build_gaussian<value_type>(target, fwhm);

    KPM_Vector<T, D> phi(2, *this);
    KPM_Vector<T, D> op_ket(1, *this);
    KPM_Vector<T, D> ket_as_kpm(1, *this);
    Eigen::Array<T, -1, 1> ket(r.Sized);
    const unsigned Orb = r.Orb;

    // map_r^O = Re(sum_{a,b} O_r(a,b) ket*_a ket_b) is a single value per
    // SITE, summed over both orbital indices -- unlike the plain diagonal
    // map, which is naturally one value per (site, orbital). results[op]
    // therefore has r.Nd (local, ghost-included, spatial-only) rows, not
    // r.Sized.
    std::vector<Eigen::Array<value_type, -1, -1>> results(
      NOperators, Eigen::Array<value_type, -1, -1>::Zero(r.Nd, 2)
    );
    Eigen::Array<value_type, -1, 1> map_op_full(r.Sized);
    Eigen::Array<value_type, -1, 1> map_op_site(r.Nd);
    Eigen::Array<value_type, -1, 1> prv(r.Nd);

    h.generate_disorder();
    for (int vec = 0; vec < vectors_; ++vec) {
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);
      ket.setZero();

      phi.Exchange_Boundaries();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }

      const value_type weight = 1.0 / (vec + 1);

      // ket_as_kpm wraps the already-propagated `ket` so multiply_orb_mtx
      // (the same on-site block-mixing routine used by custom_one/
      // custom_two and gaussian_wave_packet) can apply each operator to it
      // in one whole-lattice pass -- no per-site scalar loop.
      ket_as_kpm.v.col(0) = ket;

      for (unsigned op = 0; op < NOperators; ++op) {
        multiply_orb_mtx(operator_collection[op], &ket_as_kpm, &op_ket);
        map_op_full = factor * (ket.conjugate() * op_ket.v.col(0).array()).real();

        // Sum the b-orbital readout across the site's orbital block --
        // map_op_full is one value per (site, orbital b); the operator-
        // weighted quantity is one value per site, summed over b.
        map_op_site.setZero();
        for (unsigned b = 0; b < Orb; ++b)
          map_op_site += map_op_full.segment(b * r.Nd, r.Nd);

        prv = results[op].col(0);
        results[op].col(0) += weight * (map_op_site - results[op].col(0));
        results[op].col(1) +=
          weight * ((map_op_site - prv) * (map_op_site - results[op].col(0)) - results[op].col(1));
      }
    }

    for (unsigned op = 0; op < NOperators; ++op) {
      results[op].col(1) = results[op].col(1).sqrt() / std::sqrt(vectors_);
      store_ldos_operators(operator_labels[op], results[op]);
    }
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_ldos_operators(
  const std::string &label,
  const Eigen::Array<value_type, -1, -1> &results_
)
{
  debug_message("Entered store_ldos_operators\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
  // results_ holds one value per SITE (already summed over both orbital
  // indices, see ldos_operators above) -- not one per (site, orbital) like
  // the plain Map, so the orbital coordinate is fixed at 0 throughout: it
  // only serves to resolve the site's local/global flat index (orbital
  // isn't spatially decomposed, so this doesn't restrict which sites are
  // visited), and the output has r.Sizet/r.Orb rows instead of r.Sizet.
  const std::size_t NumGlobalSites = r.Sizet / r.Orb;
#pragma omp master
  Global.ldos_map.resize(NumGlobalSites, 2);
#pragma omp barrier
  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }
  auto body = [&](const std::array<unsigned, D> &i) {
    if constexpr (D == 2)
      local.set({i[1], i[0], 0});
    else if constexpr (D == 3)
      local.set({i[2], i[1], i[0], 0});
    r.convertCoordinates(global, local);
    Global.ldos_map.row(global.index) = results_.row(local.index);
  };
  UnitCellLoop<D>::run(idx, start, final, body);
#pragma omp barrier
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> ldos_r = Global.ldos_map.real();
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    try {
      H5::Exception::dontPrint();
      file->createGroup("/Calculation/ldos_map/Map_Operators");
    } catch (H5::Exception &e) {
      // group already exists (a previous label already created it)
    }
    write_hdf5(ldos_r, file, "/Calculation/ldos_map/Map_Operators/" + label);
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_ldos_operators\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
