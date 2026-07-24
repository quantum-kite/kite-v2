/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/


#include "Generic.hpp"
#include "BesselJ.hpp"
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

#if !(COMPILE_WAVEPACKET)
#warning "Cannot compile SimulationGaussianWavepacket.cpp. This error is not fatal, but KITE will not be able to run GaussianWavepacket(). A more recent version of gcc (8.0) is required."
#endif

template <typename T, unsigned DIM>
void Simulation<T, DIM>::calc_wavepacket(){
    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
    //

#pragma omp barrier

    // Check if the Gaussian_Wave_Packet needs to be calculated
    bool local_calculate_wavepacket;
#pragma omp master
    {
        Global.calculate_wavepacket = 0;
        H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
      try{
        int dummy_var;
        get_hdf5<int>(&dummy_var, file, (char *) "/Calculation/gaussian_wave_packet/NumDisorder");
        Global.calculate_wavepacket = 1;
      } catch(H5::Exception& e) {debug_message("Wavepacket: no need to calculate.\n");}
        file->close();  
        delete file;
    }
#pragma omp barrier
#pragma omp critical
    local_calculate_wavepacket = Global.calculate_wavepacket;
      
    // Now calculate it
#pragma omp barrier
    if(local_calculate_wavepacket){
      Gaussian_Wave_Packet();
    }
}

template <typename T, unsigned D>
void Simulation<T,D>::Gaussian_Wave_Packet(){
#if COMPILE_WAVEPACKET
#pragma omp master
      {
        std::cout << "Calculating Wavepacket.\n";
      }
#pragma omp barrier
  ComplexTraits<T> CT;
  KPM_Vector<T,D> phi (2, *this), sum_ket(1u,*this), op_ket(1u,*this);
  int NumDisorder, NumMoments, NumPoints;

  T II   = CT.assign_value(double(0), double(1));
  std::vector<double> times;

  // Observables to track during the propagation: any operator the user has
  // registered via add_orbital_index/add_orbital_coupling (spin, orbital
  // angular momentum, quadrupoles, ...), read in the order given by the
  // "Operators" label list -- not by HDF5 group-iteration order.
  std::vector<std::string> operator_labels;
  std::vector<Eigen::Matrix<std::complex<double>,-1,-1>> operator_collection;
  Eigen::Array<T,-1,-1> avg_ops, avg_ident;

  op_ket.v.col(0).setZero();
  op_ket.initiate_phases();

  float timestep;
  double width;
  Eigen::Array<T, -1, -1> avg_results;
  Eigen::Array<T, -1, -1> results(2*D, 1);
  H5::DataSet * dataset;
  H5::DataSpace * dataspace;
  hsize_t dim[2];
  Eigen::Matrix <double,-1, -1> k_vector;
  Eigen::Matrix <double ,1, D> vb;
  Eigen::Matrix <T,-1, -1>        spinor;

  //Load bra and ket
#pragma omp critical
  {
    H5::H5File * file  = new H5::H5File(name, H5F_ACC_RDONLY);
    dataset            = new H5::DataSet(file->openDataSet("/Calculation/gaussian_wave_packet/k_vector")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;

    k_vector  = Eigen::Matrix<double,-1, -1>::Zero(dim[1],dim[0]);
    spinor    = Eigen::Matrix<     T,-1, -1>::Zero(r.Orb,dim[0]);
    vb = Eigen::Matrix<     double,1, D>::Zero(D);
    get_hdf5    <int>(&NumDisorder,    file, (char *) "/Calculation/gaussian_wave_packet/NumDisorder");
    get_hdf5    <int>(&NumMoments,     file, (char *) "/Calculation/gaussian_wave_packet/NumMoments" );
    get_hdf5    <int>(&NumPoints,      file, (char *) "/Calculation/gaussian_wave_packet/NumPoints"  );
    get_hdf5  <float>(&timestep,       file, (char *) "/Calculation/gaussian_wave_packet/timestep"   );
    get_hdf5 <double>(&width,          file, (char *) "/Calculation/gaussian_wave_packet/width"      );
    get_hdf5      <T>(spinor.data(),   file, (char *) "/Calculation/gaussian_wave_packet/spinor");
    get_hdf5     <double>(vb.data(),   file, (char *) "/Calculation/gaussian_wave_packet/mean_value");
    get_hdf5 <double>(k_vector.data(), file, (char *) "/Calculation/gaussian_wave_packet/k_vector");

    // Optional: no registered operators means only Id/mean_value/Var are tracked.
    try {
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      std::string tmp = "/Calculation/gaussian_wave_packet/Operators";
      my_get_hdf5(operator_labels, my_file, tmp);
      H5::Group grp = file->openGroup("/Calculation/gaussian_wave_packet/CustomOperators/");
      for (const auto &label : operator_labels) {
        if (auto err = this->getMembers(grp, label, &operator_collection)) {
          std::cerr << "getMembers failed for " << label << "\n";
        }
      }
    } catch (H5::Exception &e) {
      debug_message("Wavepacket: no custom operators registered.\n");
    }

    file->close();  delete file;
  }
#pragma omp barrier
  const unsigned NumOperators = static_cast<unsigned>(operator_labels.size());
  avg_ops     = Eigen::Array<T,-1,-1>::Zero(NumOperators, NumPoints);
  avg_ident   = Eigen::Array<T,-1,-1>::Zero(NumPoints,1);
  avg_results = Eigen::Array<T, -1,-1>::Zero(2*D, NumPoints);

#pragma omp master
  {
    Global.avg_ops      = Eigen::Array<T,-1,-1>::Zero(NumOperators,NumPoints);
    Global.avg_ident    = Eigen::Array<T,-1,-1>::Zero(NumPoints,1);
    Global.avg_results  = Eigen::Array<T,-1,-1>::Zero(2*D,NumPoints);
  }


  NumMoments = (NumMoments/2)*2;
  Eigen::Matrix<T,-1,1> m(NumMoments);
  std::vector<double> besselj(NumMoments > 0 ? NumMoments : 0);
  if(NumMoments > 0)
    cyl_bessel_j_series(unsigned(NumMoments - 1), double(timestep), besselj.data());
  for(unsigned n = 0; n < unsigned(NumMoments); n++)
    m(n) = value_type((n == 0 ? 1 : 2 )*besselj[n]) * T(pow(-II,n));
    
  for(int id = 0; id < NumDisorder; id++)
    {
      sum_ket.set_index(0);
      sum_ket.v.setZero();
      sum_ket.build_wave_packet(k_vector, spinor, width, vb);
      h.generate_disorder();	
      sum_ket.empty_ghosts(0);
      for(unsigned t = 0; t < unsigned(NumPoints); t++)
        {
          if(t > 0)
            {
              phi.v.setZero();
              phi.set_index(0);
              phi.v.col(0) = sum_ket.v.col(0);
              phi.Exchange_Boundaries();
              phi.cheb_iteration(1); // multiply by H
              sum_ket.v.col(0) = phi.v * m.segment(0, 2);
              for(unsigned n = 2; n < unsigned(NumMoments); n += 2)
                {
                  phi.cheb_iteration(n);
                  phi.cheb_iteration(n + 1);
                  sum_ket.v.col(0) += phi.v * m.segment(n,2);
                }
            }
          sum_ket.empty_ghosts(0);

          auto x3 = sum_ket.get_point();
          avg_ident(t) += (x3 - avg_ident(t) ) /T(id + 1);

          for(unsigned op_idx = 0; op_idx < NumOperators; op_idx++)
            {
              multiply_orb_mtx(operator_collection[op_idx], &sum_ket, &op_ket);
              auto xop = sum_ket.v.adjoint() * op_ket.v.col(0);
              avg_ops(op_idx, t) += (xop(0,0) - avg_ops(op_idx, t) ) /T(id + 1);
            }

          phi.measure_wave_packet(sum_ket.v.data(), sum_ket.v.data(), results.data());
          avg_results.col(t) += (results - avg_results.col(t) ) /T(id + 1);
        }

    }

#pragma omp critical
  {
    Global.avg_ops += avg_ops;
    Global.avg_ident += avg_ident;
    Global.avg_results += avg_results;
  }
#pragma omp barrier


    
#pragma omp master
  {
    std::cout << name << std::endl;
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);

    for(unsigned op_idx = 0; op_idx < NumOperators; op_idx++)
      {
        Eigen::Array<T,-1,-1> op_series = Global.avg_ops.row(op_idx).transpose();
        std::string path = "/Calculation/gaussian_wave_packet/" + operator_labels[op_idx];
        write_hdf5(op_series, file, path);
      }
    write_hdf5(Global.avg_ident, file, (char *) "/Calculation/gaussian_wave_packet/Id");

    Eigen::Array<T,-1,-1> avg_pos(NumPoints,1);
    for(unsigned i = 0; i < D; i++)
      {
        std::string orient = "xyz";
        char name[200];
        avg_pos.col(0) = Global.avg_results.row(2*i) ;
        sprintf(name,"/Calculation/gaussian_wave_packet/mean_value%c", orient.at(i));
        write_hdf5(avg_pos, file, name);
      }

    for(unsigned i = 0; i < D; i++)
      {
        std::string orient = "xyz";
        char name[200];
        avg_pos.col(0) = Global.avg_results.row(2*i + 1) - Global.avg_results.row(2*i)*Global.avg_results.row(2*i);
        sprintf(name,"/Calculation/gaussian_wave_packet/Var%c", orient.at(i));

        write_hdf5(avg_pos, file, name);
      }

    file->close();
    delete file;
  }
#pragma omp barrier
#endif
}

#define instantiate(type,dim) template class Simulation<type,dim>;
#include "instantiate.hpp"
