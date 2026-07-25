/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include <omp.h>

#include "H5Cpp.h"
#include "../Tools/ComplexTraits.hpp"
#include "../Tools/myHDF5.hpp"

#include "../Tools/parse_input.hpp"
#include "../Tools/systemInfo.hpp"
#include "ldos.hpp"

#include "../Tools/functions.hpp"
#include "../macros.hpp"

// kite-tools' myHDF5.hpp has no vector<string> reader (unlike KITEx's, used
// by my_get_hdf5 on the simulation side) -- read the variable-length string
// list written by Python's hp.string_dtype(encoding='utf-8') directly.
static void read_string_list(std::vector<std::string> &labels, H5::H5File &file, const std::string &path){
  const H5::DataSet dataset = file.openDataSet(path);
  const H5::DataSpace dataspace = dataset.getSpace();
  const H5::StrType str_type = dataset.getStrType();
  hsize_t dim[1];
  dataspace.getSimpleExtentDims(dim, NULL);
  labels.resize(dim[0]);
  std::vector<char*> buffer(dim[0], nullptr);
  dataset.read(buffer.data(), str_type);
  for(hsize_t i = 0; i < dim[0]; i++){
    labels[i] = buffer[i];
    free(buffer[i]);
  }
}

template <typename T, unsigned DIM>
ldos<T, DIM>::ldos(system_info<T, DIM>& sysinfo, shell_input & vari){
    // Class constructor
    
    systemInfo  = &sysinfo;              // retrieve the information about the Hamiltonian
    variables   = vari;                   // retrieve the shell input
    dirName     = "/Calculation/ldos/";     // location of the information about the conductivity
    
    isRequired = is_required() && variables.lDOS_is_required;         // was the local density of states requested?
    isPossible = false;                 // do we have all we need to calculate the density of states?
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();

      if(isPossible){
          printLDOS();                  // Print all the parameters used
          calculate();
          calculate_operators();
      } else {
        std::cout << "ERROR. The LDOS was requested but the data "
            "needed for its computation was not found in the input .h5 file. "
            "Make sure KITEx has processed the file first. Exiting.";
        exit(1);
      }
    }
}

template <typename T, unsigned DIM>
void ldos<T, DIM>::printLDOS(){
  double scale = systemInfo->energy_scale;
    std::cout << "The local density of states will be calculated with the following parameters:\n"
        "   Number of energies: " << NumEnergies << "\n"
        "   Number of positions: " << NumPositions << "\n"
        "   Filename: " << filename  << "X.dat" << ((default_filename)?" (default)":"") << "\n"
        "   Kernel: "               << kernel           << ((default_kernel)?           " (default)":"") << "\n";
    if(kernel == "green"){
        std::cout << "   Kernel parameter: "     << kernel_parameter*scale << ((default_kernel_parameter)? " (default)":"") << "\n";
    }
}

template <typename T, unsigned DIM>
bool ldos<T, DIM>::is_required(){
    // check whether the local density of states was asked for
    // if this quantity exists, so should all the others.

    name = systemInfo->filename;
	H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    bool result = false;
    try{
        H5::Exception::dontPrint();
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());									
        result = true;
    } catch(H5::Exception& e){}
  

    file.close();

    return result;
}
	
template <typename T, unsigned DIM>
void ldos<T, DIM>::override_parameters(){
    if(variables.lDOS_Name != ""){
        filename         = variables.lDOS_Name;
        default_filename = false;
    }

    if(variables.lDOS_NumMoments != -1){
        NumMoments         = variables.lDOS_NumMoments;
        default_NumMoments = false;
        if(variables.lDOS_NumMoments > MaxMoments){
          std::cout << "lDOS: The number of Chebyshev moments specified"
            " cannot be larger than the number of moments calculated by KITEx."
            " Please specify a smaller number. Exiting.\n";
          exit(1);
        }
    }

    //std::cout << "variables kernel:" << variables.lDOS_kernel << "\n";
    if(variables.lDOS_kernel != ""){
      //std::cout << "entered if \n";
        kernel         = variables.lDOS_kernel;
        default_kernel = false;
    }
    //std::cout << "kernel: " << kernel << "\n";

    if(kernel == "green"){
      if(variables.lDOS_kernel_parameter != -8888.8){
        kernel_parameter = variables.lDOS_kernel_parameter/systemInfo->energy_scale;
        default_kernel_parameter = false;
      }
    }
}

template <typename T, unsigned DIM>
void ldos<T, DIM>::set_default_parameters(){
    filename = "ldos";
    default_filename = true;
    MaxMoments = -1;

    // kernel options
    kernel = "jackson";
    default_kernel = true;
    default_kernel_parameter = true;

}


template <typename T, unsigned DIM>
bool ldos<T, DIM>::fetch_parameters(){
  debug_message("Entered ldos::fetch_parameters.\n");
  //This function reads all the data from the hdf5 file that's needed to 
  //calculate the LDoS
  
  // Check if the data for the ldos exists
  if(!isRequired){
    std::cout << "Data for LDoS does not exist. Exiting.\n";
    exit(1);
  }
  
  H5::DataSet * dataset;
  H5::DataSpace * dataspace;
  hsize_t dim[2];
  H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
  
  dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Orbitals")  );
  dataspace          = new H5::DataSpace(dataset->getSpace());
  dataspace -> getSimpleExtentDims(dim, NULL);
  dataspace->close(); delete dataspace;
  dataset->close();   delete dataset;
  NumPositions = dim[0];
  
  dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Energy")  );
  dataspace          = new H5::DataSpace(dataset->getSpace());
  dataspace -> getSimpleExtentDims(dim, NULL);
  dataspace->close(); delete dataspace;
  dataset->close();   delete dataset;
  NumEnergies = dim[0];
  
  ldos_Orbitals = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,1);
  ldos_Positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,1);
  energies = Eigen::Matrix<float, -1, -1>::Zero(NumEnergies,1);
  
  //Fetch the relevant parameters from the hdf file
  get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());	
  get_hdf5(ldos_Orbitals.data(), &file, (char*)"/Calculation/ldos/Orbitals");
  get_hdf5(ldos_Positions.data(), &file, (char*)"/Calculation/ldos/FixPosition");
  get_hdf5(energies.data(), &file, (char*)"/Calculation/ldos/Energy");
  
  if(DIM == 2){
    global_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,3);
    for(long i = 0; i < NumPositions; i++){
      int Lx = systemInfo->size[0];
      global_positions(i,0) = ldos_Positions(i)%Lx;
      global_positions(i,1) = ldos_Positions(i)/Lx;
      global_positions(i,2) = ldos_Orbitals(i);
    }
  } else if(DIM ==3){
    global_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,4);
    for(long i = 0; i < NumPositions; i++){
      int Lx = systemInfo->size[0];
      int Ly = systemInfo->size[1];
      
      global_positions(i,0) = ldos_Positions(i)%(Lx);
      global_positions(i,1) = (ldos_Positions(i)%(Lx*Ly))/Lx;
      global_positions(i,2) = ldos_Positions(i)/(Lx*Ly);
      global_positions(i,3) = ldos_Orbitals(i);
    }
  }
  
  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo->isComplex;
  
  bool result = false;
  // Retrieve the lmu Matrix
  std::string MatrixName = dirName + "lMU";
  try{
    debug_message("Filling the lMU matrix.\n");
    lMU = Eigen::Matrix<std::complex<T>,-1,-1>::Zero(MaxMoments, NumPositions);
    
    if(complex)
      get_hdf5(lMU.data(), &file, (char*)MatrixName.c_str());
    if(!complex){
      Eigen::Matrix<T,-1,-1> lMUReal; 
      lMUReal = Eigen::Matrix<T,-1,-1>::Zero(MaxMoments, NumPositions); 
      get_hdf5(lMUReal.data(), &file, (char*)MatrixName.c_str()); 
      
      lMU = lMUReal.template cast<std::complex<T>>();
    }				
    
    result = true;
  } catch(H5::Exception& e) {debug_message("lDOS: There is no lMU matrix.\n");}


  NumMoments = MaxMoments;

  // Optional: operator-weighted moments (Tr[O*Im G(r,r,E)]), written by
  // Src/Simulation/Custom/SimulationLDOSOperator.cpp only when
  // calculation.ldos(..., operators=[...]) was used. Absent otherwise.
  try{
    read_string_list(operator_labels, file, dirName + "Operators");

    H5::DataSet * opdataset;
    H5::DataSpace * opdataspace;
    hsize_t opdim[1];
    opdataset   = new H5::DataSet(file.openDataSet(dirName + "OperatorPositions"));
    opdataspace = new H5::DataSpace(opdataset->getSpace());
    opdataspace -> getSimpleExtentDims(opdim, NULL);
    opdataspace->close(); delete opdataspace;
    opdataset->close();   delete opdataset;
    NumOperatorPositions = opdim[0];

    op_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumOperatorPositions,1);
    get_hdf5(op_positions.data(), &file, (char*)(dirName+"OperatorPositions").c_str());

    if(DIM == 2){
      global_op_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumOperatorPositions,2);
      for(unsigned i = 0; i < NumOperatorPositions; i++){
        int Lx = systemInfo->size[0];
        global_op_positions(i,0) = op_positions(i)%Lx;
        global_op_positions(i,1) = op_positions(i)/Lx;
      }
    } else if(DIM == 3){
      global_op_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumOperatorPositions,3);
      for(unsigned i = 0; i < NumOperatorPositions; i++){
        int Lx = systemInfo->size[0];
        int Ly = systemInfo->size[1];
        global_op_positions(i,0) = op_positions(i)%(Lx);
        global_op_positions(i,1) = (op_positions(i)%(Lx*Ly))/Lx;
        global_op_positions(i,2) = op_positions(i)/(Lx*Ly);
      }
    }

    for(const auto &label : operator_labels){
      Eigen::Matrix<std::complex<T>,-1,-1> mat =
        Eigen::Matrix<std::complex<T>,-1,-1>::Zero(MaxMoments, NumOperatorPositions);
      std::string opMatrixName = dirName + "lMU_Operators/" + label;
      if(complex){
        get_hdf5(mat.data(), &file, (char*)opMatrixName.c_str());
      } else {
        Eigen::Matrix<T,-1,-1> matReal =
          Eigen::Matrix<T,-1,-1>::Zero(MaxMoments, NumOperatorPositions);
        get_hdf5(matReal.data(), &file, (char*)opMatrixName.c_str());
        mat = matReal.template cast<std::complex<T>>();
      }
      lMU_Operators.push_back(mat);
    }
  } catch(H5::Exception& e) {
    debug_message("lDOS: no operator-weighted moments requested.\n");
    operator_labels.clear();
    lMU_Operators.clear();
  }

  file.close();
  debug_message("Left lDOS::fetch_parameters.\n");
  return result;
}


template <typename U, unsigned DIM>
void ldos<U, DIM>::calculate(){
  
  Eigen::Matrix<std::complex<U>, -1, -1> LDOS;
  LDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);
  
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> OrderedMU;
  OrderedMU = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(NumEnergies, NumPositions);
  OrderedMU = lMU;
  
  omp_set_num_threads(systemInfo->NumThreads);
  //omp_set_num_threads(1);
#pragma omp parallel 
  {
#pragma omp critical
    {
      int localN = NumMoments/systemInfo->NumThreads;
      //int localN = NumMoments;
      int thread_id = omp_get_thread_num();
      //std::cout << "thread_id: " << thread_id << "\n";
      //std::cout << "NumMoments: " << localN << "\n";
      //thread_id = 0;
      long offset = thread_id*localN*NumPositions;
      Eigen::Map<Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>> locallMU(OrderedMU.data() + offset, localN, NumPositions);
      
      //std::cout << "locallMU: \n" << locallMU << "\n";
      
      Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
      GammaE = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, localN);
      
      U factor;
      
      if(kernel == "jackson"){
	for(int i = 0; i < NumEnergies; i++){
	  for(int m = 0; m < localN; m++){
	    factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
	    GammaE(i,m) += delta(m + thread_id*localN,energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
	  }
	}
      }
      
      
      if(kernel == "green"){
	std::complex<U> c_energy;
	for(int i = 0; i < NumEnergies; i++){
	  c_energy = std::complex<U>(energies(i), kernel_parameter);
	  for(int m = 0; m < localN; m++){
	    factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
	    GammaE(i,m) += -factor*green<std::complex<U>>(m, 1, c_energy).imag();
	  }
	}
      }
      
      //std::cout << "GammaE: \n" << GammaE << "\n";
      //for(int m = 0; m < ; m++){
      //factor = 1.0/(1.0 + U(m==0));
      //for(int i = 0; i < NumEnergies; i++){
      //GammaE(i,m) = delta(m + thread_id*localN, energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
      //}
      //}
      
      Eigen::Matrix<std::complex<U>, -1, -1> localLDOS;
      localLDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);
      localLDOS = GammaE*locallMU;
      //#pragma omp critical
      LDOS += localLDOS;
    }
  }
  
  // Save the density of states to a file
  U mult = 1.0/systemInfo->energy_scale;
  std::ofstream myfile;
  double scale = systemInfo->energy_scale;
  double shift = systemInfo->energy_shift;
  for(int i=0; i < NumEnergies; i++){
    myfile.open(filename + std::to_string(energies(i)*scale + shift) + ".dat");
    if(DIM == 2){
      for(unsigned pos = 0; pos < NumPositions; pos++){
	int x, y, orb;
	x = global_positions(pos,0);
	y = global_positions(pos,1);
	orb = global_positions(pos,2);
	myfile  << x << " " << y << " " << orb << " " << LDOS(i,pos).real()*mult << "\n";
      };
    } else if(DIM == 3){
      for(unsigned pos = 0; pos < NumPositions; pos++){
	int x, y, z, orb;
	x = global_positions(pos,0);
	y = global_positions(pos,1);
	z = global_positions(pos,2);
	orb = global_positions(pos,3);
	myfile  << x << " " << y << " " << z << " " << orb << " " << LDOS(i,pos).real()*mult << "\n";
      };
    }
    myfile.close();
  }
}


template <typename U, unsigned DIM>
void ldos<U, DIM>::calculate_operators(){
  // Reconstructs Tr[O*Im G(r,r,E)] for each registered operator, reusing the
  // exact same Jackson/Green-kernel Chebyshev reconstruction as calculate()
  // above -- only the moment matrix and position list differ, since the
  // reconstruction stage is operator-agnostic.
  if(operator_labels.empty()) return;

  for(std::size_t op_idx = 0; op_idx < operator_labels.size(); op_idx++){
    const std::string &label = operator_labels[op_idx];

    Eigen::Matrix<std::complex<U>, -1, -1> LDOS;
    LDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumOperatorPositions);

    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> OrderedMU;
    OrderedMU = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(NumEnergies, NumOperatorPositions);
    OrderedMU = lMU_Operators[op_idx];

    omp_set_num_threads(systemInfo->NumThreads);
#pragma omp parallel
    {
#pragma omp critical
      {
        int localN = NumMoments/systemInfo->NumThreads;
        int thread_id = omp_get_thread_num();
        long offset = thread_id*localN*NumOperatorPositions;
        Eigen::Map<Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>> locallMU(OrderedMU.data() + offset, localN, NumOperatorPositions);

        Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
        GammaE = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, localN);

        U factor;

        if(kernel == "jackson"){
          for(int i = 0; i < NumEnergies; i++){
            for(int m = 0; m < localN; m++){
              factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
              GammaE(i,m) += delta(m + thread_id*localN,energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
            }
          }
        }

        if(kernel == "green"){
          std::complex<U> c_energy;
          for(int i = 0; i < NumEnergies; i++){
            c_energy = std::complex<U>(energies(i), kernel_parameter);
            for(int m = 0; m < localN; m++){
              factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
              GammaE(i,m) += -factor*green<std::complex<U>>(m, 1, c_energy).imag();
            }
          }
        }

        Eigen::Matrix<std::complex<U>, -1, -1> localLDOS;
        localLDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumOperatorPositions);
        localLDOS = GammaE*locallMU;
        LDOS += localLDOS;
      }
    }

    U mult = 1.0/systemInfo->energy_scale;
    std::ofstream myfile;
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;
    for(int i=0; i < NumEnergies; i++){
      myfile.open(filename + "_" + label + "_" + std::to_string(energies(i)*scale + shift) + ".dat");
      if(DIM == 2){
        for(unsigned pos = 0; pos < NumOperatorPositions; pos++){
          int x, y;
          x = global_op_positions(pos,0);
          y = global_op_positions(pos,1);
          myfile  << x << " " << y << " " << LDOS(i,pos).real()*mult << "\n";
        };
      } else if(DIM == 3){
        for(unsigned pos = 0; pos < NumOperatorPositions; pos++){
          int x, y, z;
          x = global_op_positions(pos,0);
          y = global_op_positions(pos,1);
          z = global_op_positions(pos,2);
          myfile  << x << " " << y << " " << z << " " << LDOS(i,pos).real()*mult << "\n";
        };
      }
      myfile.close();
    }
  }
}


// Instantiations
template class ldos<float, 1u>;
template class ldos<float, 2u>;
template class ldos<float, 3u>;

template class ldos<double, 1u>;
template class ldos<double, 2u>;
template class ldos<double, 3u>;

template class ldos<long double, 1u>;
template class ldos<long double, 2u>;
template class ldos<long double, 3u>;
