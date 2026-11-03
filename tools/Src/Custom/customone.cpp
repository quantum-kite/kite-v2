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
#include "customone.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
customone<T, DIM>::customone(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    customone_finished = false;
    MaxMoments = -1;
    isPossible = false;
    isRequired = is_required() && variables.CustomOne_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        energies = Eigen::Matrix<T, -1, 1>::LinSpaced(NEnergies, Emin, Emax);
        if(isPossible){
            printCustomOne();
            calculate();
        } else {
            std::cout << "ERROR. custom_one() post-processing was requested but the data "
                "needed for its computation (/Calculation/CustomOne/*) was not found in "
                "the input .h5 file. Make sure KITEx has processed the file first, and "
                "that calculation.custom_one(...) was actually called. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool customone<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomOne/";
    bool result = false;
    try{
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void customone<T, DIM>::set_default_parameters(){
    NEnergies = 1024;
    default_NEnergies = true;

    filename = "custom_one.dat";
    default_filename = true;

    Emax = 0.99;
    Emin = -0.99;
    default_Emin = true;
    default_Emax = true;

    kernel = "jackson";
    default_kernel = true;
    default_kernel_parameter = true;
}

template <typename T, unsigned DIM>
void customone<T, DIM>::override_parameters(){
    // Overrides the current parameters with the ones from the shell input.
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    if(variables.CustomOne_NumEnergies != -1){
        NEnergies         = variables.CustomOne_NumEnergies;
        default_NEnergies = false;
    }

    if(variables.CustomOne_Emin != 8888.8){
        Emin = (variables.CustomOne_Emin - shift)/scale;
        default_NEnergies = false;
    }

    if(variables.CustomOne_Emax != -8888.8){
        Emax = (variables.CustomOne_Emax - shift)/scale;
        default_NEnergies = false;
    }

    if(variables.CustomOne_Name != ""){
        filename         = variables.CustomOne_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void customone<T, DIM>::printCustomOne(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;
    std::string energy_range = "[" + std::to_string(Emin*scale + shift) + ", " + std::to_string(Emax*scale + shift) + "]";

    std::cout << "The custom_one() spectral density will be calculated with these parameters: (eV)\n"
        "   Energy range: "        << energy_range << ((default_Emin && default_Emax)? " (default)":"") << "\n"
        "   Number of energies: "  << NEnergies     << ((default_NEnergies)?           " (default)":"") << "\n"
        "   Filename: "            << filename      << ((default_filename)?            " (default)":"") << "\n"
        "   Number of moments: "   << NumMoments     << "\n"
        "   Number of velocity operators in the vertex: " << NumVelocities <<
        ((NumVelocities % 2 == 0)? " (even -> Tr[T_n(H)*A] is real)"
                                  : " (odd -> Tr[T_n(H)*A] is imaginary)") << "\n"
        "   Kernel: jackson\n";
}

template <typename T, unsigned DIM>
bool customone<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomOne/";

    get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());
    NumMoments = MaxMoments;

    get_hdf5(&NumVelocities, &file, (char*)(dirName+"NumVelocities").c_str());

    // Gamma is always complex128 (custom_one requires is_complex=True) -- unlike
    // dos<T,DIM>'s MU, which can be real or complex depending on the Hamiltonian.
    bool result = false;
    std::string MatrixName = dirName + "Gamma";
    try{
        Gamma = Eigen::Array<std::complex<T>,-1,-1>::Zero(1, MaxMoments);
        get_hdf5(Gamma.data(), &file, (char*)MatrixName.c_str());
        result = true;
    } catch(H5::Exception& e) {debug_message("CustomOne: There is no Gamma matrix.\n");}

    file.close();
    return result;
}

template <typename U, unsigned DIM>
void customone<U, DIM>::calculate(){
    GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

    U scale = systemInfo->energy_scale;
    U mult = 1.0/scale;
    U factor;
    U shift = systemInfo->energy_shift;

    // store_custom_one (Src/Simulation/Custom/SimulationRankOne.cpp) already keeps
    // whichever component of the raw stochastic trace is physically meaningful --
    // Re(gamma) for an even NumVelocities, or i*Im(gamma) (still complex-valued,
    // purely imaginary) for an odd one. Rotate the odd case back onto the real
    // axis (multiply by -i, i.e. divide by i) before feeding it through the exact
    // same real-valued Jackson-kernel reconstruction dos<T,DIM> uses -- this is the
    // ONLY place the "i" bookkeeping shows up in post-processing; everything else
    // below is identical to dos.cpp.
    const std::complex<U> rotation = (NumVelocities % 2 == 0)
        ? std::complex<U>(1.0, 0.0)
        : std::complex<U>(0.0, -1.0);

    for(int i = 0; i < NEnergies; i++){
        for(int m = 0; m < NumMoments; m++){
            factor = 1.0/(1.0 + U(m==0));
            GammaE(i) += (rotation*Gamma(m))*delta(m,energies(i))*kernel_jackson<U>(m, NumMoments)*factor*mult;
        }
    }

    std::ofstream myfile;
    myfile.open(filename);
    for(int i=0; i < NEnergies; i++){
        myfile << energies(i)*scale + shift << " " << GammaE.real()(i) << "\n";
    }
    myfile.close();
    customone_finished = true;
}

template class customone<float, 1u>;
template class customone<float, 2u>;
template class customone<float, 3u>;

template class customone<double, 1u>;
template class customone<double, 2u>;
template class customone<double, 3u>;

template class customone<long double, 1u>;
template class customone<long double, 2u>;
template class customone<long double, 3u>;
