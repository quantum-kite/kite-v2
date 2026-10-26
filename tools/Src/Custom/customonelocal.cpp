/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include "H5Cpp.h"
#include "../Tools/ComplexTraits.hpp"
#include "../Tools/myHDF5.hpp"

#include "../Tools/parse_input.hpp"
#include "../Tools/systemInfo.hpp"
#include "customonelocal.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
customonelocal<T, DIM>::customonelocal(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    customonelocal_finished = false;
    MaxMoments = -1;
    NumPositions = 0;
    isPossible = false;
    isRequired = is_required() && variables.CustomOneLocal_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        energies = Eigen::Matrix<T, -1, 1>::LinSpaced(NEnergies, Emin, Emax);
        if(isPossible){
            printCustomOneLocal();
            calculate();
        } else {
            std::cout << "ERROR. custom_one_local() post-processing was requested but the "
                "data needed for its computation (/Calculation/CustomOneLocal/*) was not "
                "found in the input .h5 file. Make sure KITEx has processed the file first, "
                "and that calculation.custom_one_local(...) was actually called. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool customonelocal<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomOneLocal/";
    bool result = false;
    try{
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void customonelocal<T, DIM>::set_default_parameters(){
    NEnergies = 1024;
    default_NEnergies = true;

    filename = "custom_one_local.dat";
    default_filename = true;

    Emax = 0.99;
    Emin = -0.99;
    default_Emin = true;
    default_Emax = true;
}

template <typename T, unsigned DIM>
void customonelocal<T, DIM>::override_parameters(){
    // Overrides the current parameters with the ones from the shell input.
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    if(variables.CustomOneLocal_NumEnergies != -1){
        NEnergies         = variables.CustomOneLocal_NumEnergies;
        default_NEnergies = false;
    }

    if(variables.CustomOneLocal_Emin != 8888.8){
        Emin = (variables.CustomOneLocal_Emin - shift)/scale;
        default_Emin = false;
    }

    if(variables.CustomOneLocal_Emax != -8888.8){
        Emax = (variables.CustomOneLocal_Emax - shift)/scale;
        default_Emax = false;
    }

    if(variables.CustomOneLocal_Name != ""){
        filename         = variables.CustomOneLocal_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void customonelocal<T, DIM>::printCustomOneLocal(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;
    std::string energy_range = "[" + std::to_string(Emin*scale + shift) + ", " + std::to_string(Emax*scale + shift) + "]";

    std::cout << "The custom_one_local() spectral density will be calculated with these parameters: (eV)\n"
        "   Energy range: "        << energy_range  << ((default_Emin && default_Emax)? " (default)":"") << "\n"
        "   Number of energies: "  << NEnergies     << ((default_NEnergies)?            " (default)":"") << "\n"
        "   Filename: "            << filename      << ((default_filename)?             " (default)":"") << "\n"
        "   Number of moments: "   << NumMoments    << "\n"
        "   Number of probed orbitals: " << NumPositions << "\n"
        "   Kernel: jackson\n";
}

template <typename T, unsigned DIM>
bool customonelocal<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomOneLocal/";

    // The number of probed orbitals is the length of Orbitals; it is not stored
    // anywhere as a scalar.
    H5::DataSet * dataset;
    H5::DataSpace * dataspace;
    hsize_t dim[2];
    dataset   = new H5::DataSet(file.openDataSet(dirName + "Orbitals"));
    dataspace = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;
    NumPositions = dim[0];

    get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());
    NumMoments = MaxMoments;

    Orbitals  = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions, 1);
    Positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions, 1);
    get_hdf5(Orbitals.data(),  &file, (char*)(dirName+"Orbitals").c_str());
    get_hdf5(Positions.data(), &file, (char*)(dirName+"FixPosition").c_str());

    // Gamma is stored as (NumPositions, NumMoments) row-major, which is the same
    // memory as an Eigen (NumMoments, NumPositions) column-major array -- the same
    // convention ldos<T,DIM> uses for lMU, and verified against it: with an
    // identity vertex the two arrays agree exactly.
    bool result = false;
    std::string MatrixName = dirName + "Gamma";
    try{
        Gamma = Eigen::Array<std::complex<T>,-1,-1>::Zero(MaxMoments, NumPositions);
        if(systemInfo->isComplex){
            get_hdf5(Gamma.data(), &file, (char*)MatrixName.c_str());
        } else {
            Eigen::Array<T,-1,-1> GammaReal =
                Eigen::Array<T,-1,-1>::Zero(MaxMoments, NumPositions);
            get_hdf5(GammaReal.data(), &file, (char*)MatrixName.c_str());
            Gamma = GammaReal.template cast<std::complex<T>>();
        }
        result = true;
    } catch(H5::Exception& e) {debug_message("CustomOneLocal: There is no Gamma matrix.\n");}

    file.close();
    return result;
}

template <typename U, unsigned DIM>
void customonelocal<U, DIM>::calculate(){
    GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, NumPositions);

    U scale = systemInfo->energy_scale;
    U mult = 1.0/scale;
    U factor;
    U shift = systemInfo->energy_shift;

    // Identical to customone<T,DIM>::calculate(), minus the velocity-parity
    // rotation (no NumVelocities is stored for the Local variants) and minus any
    // 1/N (the seeding is deterministic, not a stochastic trace). The inner sum
    // is the standard Jackson-damped Chebyshev reconstruction of a spectral
    // density; delta() carries the 2/(pi*sqrt(1-Etilde^2)) weight.
    for(long j = 0; j < NumPositions; j++){
        for(int i = 0; i < NEnergies; i++){
            for(int m = 0; m < NumMoments; m++){
                factor = 1.0/(1.0 + U(m==0));
                GammaE(i,j) += Gamma(m,j)*delta(m,energies(i))*kernel_jackson<U>(m, NumMoments)*factor*mult;
            }
        }
    }

    std::ofstream myfile;
    myfile.open(filename);
    myfile << "# Local operator-weighted spectral density at each probed orbital,\n"
              "#   f_A(j,E) = (1/2) <j| { delta(E-H), A } |j>\n"
              "# reconstructed from Gamma_n = <j|T_n(Htilde) A|j> with the Jackson kernel.\n"
              "# No 1/N: the seeding is a single basis ket, not a stochastic trace.\n"
              "# Units: [A] per eV per orbital. For a vertex diagonal in the orbital basis\n"
              "# this is A_jj times the orbital-resolved LDoS.\n"
              "# col 1: E (eV)\n";
    for(long j = 0; j < NumPositions; j++)
        myfile << "# col " << j+2 << ": cell " << Positions(j) << ", orbital " << Orbitals(j) << "\n";
    for(int i=0; i < NEnergies; i++){
        myfile << energies(i)*scale + shift;
        for(long j = 0; j < NumPositions; j++)
            myfile << " " << GammaE.real()(i,j);
        myfile << "\n";
    }
    myfile.close();
    customonelocal_finished = true;
}

template class customonelocal<float, 1u>;
template class customonelocal<float, 2u>;
template class customonelocal<float, 3u>;

template class customonelocal<double, 1u>;
template class customonelocal<double, 2u>;
template class customonelocal<double, 3u>;

template class customonelocal<long double, 1u>;
template class customonelocal<long double, 2u>;
template class customonelocal<long double, 3u>;
