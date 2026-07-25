/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Operator-weighted generalization of the exact/deterministic LDOS method
// (SimulationLMU.cpp): instead of the scalar diagonal moment
// mu_n(r) = <r|T_n(H)|r>, computes, for each user-registered on-site
// operator O_r, mu_n^O(r) = sum_{a,b} (O_r)_{b,a} * <r,b|T_n(H)|r,a>.
// Kept in its own translation unit (mirrors custom_one/custom_two) so the
// plain LMU/calc_LDOS path is untouched.

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
void Simulation<T, D>::store_LMU_Operators(
  const std::string &label,
  Eigen::Array<T, -1, -1> *gamma
)
{
  debug_message("Entered store_LMU_Operators\n");

  long int nMoments   = gamma->rows();
  long int nPositions = gamma->cols();

#pragma omp master
  Global.general_gamma = Eigen::Array<T, -1, -1>::Zero(nMoments, nPositions);
#pragma omp barrier
#pragma omp critical
  Global.general_gamma += *gamma;
#pragma omp barrier

#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    try {
      H5::Exception::dontPrint();
      file->createGroup("/Calculation/ldos/lMU_Operators");
    } catch (H5::Exception &e) {
      // group already exists (a previous label already created it)
    }
    write_hdf5(Global.general_gamma, file, "/Calculation/ldos/lMU_Operators/" + label);
    file->close();
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_LMU_Operators\n");
}

template <typename T, unsigned D>
void Simulation<T, D>::LMU_Operators(
  int NDisorder,
  int NMoments,
  const Eigen::Array<unsigned long, -1, 1> &positions,
  const std::vector<std::string> &operator_labels,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operator_collection
)
{
  debug_message("Entered Simulation::LMU_Operators\n");

  typedef typename extract_value_type<T>::value_type value_type;
  const int NPositions = positions.size();
  const unsigned NOperators = static_cast<unsigned>(operator_labels.size());
  const unsigned Orb = r.Orb;

  unsigned long totalSites;
  if constexpr (D == 2)
    totalSites = static_cast<unsigned long>(r.Lt[0]) * r.Lt[1];
  else if constexpr (D == 3)
    totalSites = static_cast<unsigned long>(r.Lt[0]) * r.Lt[1] * r.Lt[2];
  else
    totalSites = static_cast<unsigned long>(r.Lt[0]);

  KPM_Vector<T, D> kpm0(1, *this); // seed vector, one orbital at a time

  std::vector<Eigen::Array<T, -1, -1>> gamma(
    NOperators, Eigen::Array<T, -1, -1>::Zero(NMoments, NPositions)
  );
  Eigen::Array<long, -1, 1> average = Eigen::Array<long, -1, 1>::Zero(NPositions, 1);
  Eigen::Array<T, -1, -1> contrib(NOperators, NMoments);

  for (int disorder = 0; disorder < NDisorder; disorder++) {
    h.generate_disorder();
    h.generate_twists();

    for (int pos_index = 0; pos_index < NPositions; pos_index++) {
      const unsigned long site_pos = positions(pos_index);
      std::size_t base_local_index;
      const bool correct_thread = kpm0.locate_site(site_pos, base_local_index);

      contrib.setZero();

      for (unsigned a = 0; a < Orb; a++) {
        const unsigned long pos_a = site_pos + static_cast<unsigned long>(a) * totalSites;
        kpm0.build_site(pos_a);

        // A fresh KPM_Vector per seed, not a reused/reset one: reusing a
        // single kpm1 object across many (position, orbital) seeds within
        // this loop left some internal state uncleared by
        // set_index()+v.col(0)=seed+Exchange_Boundaries() alone, causing
        // the propagated norm to grow compounding across seeds (empirically
        // confirmed: reused-object norms grew ~1.8x per seed, unboundedly;
        // a fresh object per seed stays correctly bounded). initiate_phases()
        // is safe to call again here -- it only recomputes the boundary
        // twist factors from h.BoundTwist, which is fixed for this disorder
        // realization, not re-randomized.
        KPM_Vector<T, D> kpm1(2, *this);
        kpm1.initiate_phases();
        kpm1.set_index(0);
        kpm1.v.col(0) = kpm0.v.col(0);
        kpm1.Exchange_Boundaries();
        kpm0.empty_ghosts(0);

        for (int n = 0; n < NMoments; n += 2) {
          kpm1.cheb_iteration(n);
          kpm1.cheb_iteration(n + 1);

          if (correct_thread) {
            for (unsigned b = 0; b < Orb; b++) {
              const std::size_t local_b = base_local_index + b * r.Nd;
              const T v0 = kpm1.v(local_b, 0);
              const T v1 = kpm1.v(local_b, 1);
              for (unsigned op = 0; op < NOperators; op++) {
                const auto &orb_mtx = operator_collection[op];
                // mu_n^O(r) = Tr[O * M_n(r)] = sum_{a,b} O_{ab} * M_n(b,a),
                // where M_n(b,a) = <r,b|T_n(H)|r,a> is exactly what v0/v1
                // hold here (seeded at orbital a, read out at orbital b).
                // The coefficient is therefore O_{ab} = orb_mtx(a,b), NOT
                // orb_mtx(b,a) -- using the transposed index is invisible
                // for symmetric/diagonal O (Sz, real projectors) but flips
                // the sign for antisymmetric imaginary O (e.g. Sy).
                T coeff;
                if constexpr (std::is_same<T, value_type>::value)
                  coeff = orb_mtx(a, b).real();
                else
                  coeff = static_cast<T>(orb_mtx(a, b));
                contrib(op, n)     += coeff * v0;
                contrib(op, n + 1) += coeff * v1;
              }
            }
          }
        }
      }

      for (unsigned op = 0; op < NOperators; op++)
        for (int n = 0; n < NMoments; n++)
          gamma[op](n, pos_index) +=
            (contrib(op, n) - gamma[op](n, pos_index)) / value_type(average(pos_index) + 1);
      average(pos_index)++;
    }
  }

  for (unsigned op = 0; op < NOperators; op++)
    store_LMU_Operators(operator_labels[op], &gamma[op]);

  debug_message("Left Simulation::LMU_Operators\n");
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_LDOS_operators()
{
  debug_message("Entered Simulation::calc_LDOS_operators\n");
#pragma omp barrier

  std::vector<std::string> operator_labels;
  std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> operator_collection;
  unsigned NMoments = 0, NDisorder = 0;
  Eigen::Array<unsigned long, -1, 1> op_positions;
  bool has_operators = false;

#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    try {
      std::string tmp = "/Calculation/ldos/Operators";
      my_get_hdf5(operator_labels, *file, tmp);

      get_hdf5<unsigned>(&NMoments, file, (char *) "/Calculation/ldos/NumMoments");
      get_hdf5<unsigned>(&NDisorder, file, (char *) "/Calculation/ldos/NumDisorder");

      H5::DataSet *dataset;
      H5::DataSpace *dataspace;
      hsize_t dim[1];
      dataset   = new H5::DataSet(file->openDataSet("/Calculation/ldos/OperatorPositions"));
      dataspace = new H5::DataSpace(dataset->getSpace());
      dataspace->getSimpleExtentDims(dim, NULL);
      dataspace->close(); delete dataspace;
      dataset->close();   delete dataset;

      op_positions = Eigen::Array<unsigned long, -1, 1>::Zero(dim[0], 1);
      get_hdf5<unsigned long>(op_positions.data(), file, (char *) "/Calculation/ldos/OperatorPositions");

      H5::Group grp = file->openGroup("/Calculation/ldos/CustomOperators/");
      for (const auto &label : operator_labels) {
        if (auto err = this->getMembers(grp, label, &operator_collection)) {
          std::cerr << "getMembers failed for " << label << "\n";
        }
      }
      has_operators = !operator_labels.empty();
    } catch (H5::Exception &e) {
      debug_message("ldos: no operator-weighted LDOS requested.\n");
    }
    file->close();
    delete file;
    Global.calculate_ldos_operators = has_operators;
  }
#pragma omp barrier

  bool local_has_operators;
#pragma omp critical
  local_has_operators = Global.calculate_ldos_operators;
#pragma omp barrier

  if (local_has_operators) {
#pragma omp master
    std::cout << "Calculating operator-weighted LDoS.\n";
#pragma omp barrier
    LMU_Operators(NDisorder, NMoments, op_positions, operator_labels, operator_collection);
  }
  debug_message("Left Simulation::calc_LDOS_operators\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
