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
#include "customtwo.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
customtwo<T, DIM>::customtwo(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    isPossible = false;
    isRequired = is_required() && variables.CustomTwo_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        if(isPossible){
            printCustomTwo();
            calculate();
        } else {
            std::cout << "ERROR. custom_two() post-processing was requested but the data "
                "needed for its computation (/Calculation/CustomTwo/*) was not found in "
                "the input .h5 file, or the two vertices did not have equal NumMoments "
                "(this reconstruction requires equal moment counts on both vertices, same "
                "as store_custom_two's Hermitization step). Make sure KITEx has processed "
                "the file first, and that calculation.custom_two(...) was actually called "
                "with stream_=[A,B] using the same number of moments for A and B. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool customtwo<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomTwo/";
    bool result = false;
    try{
        int dummy;
        get_hdf5(&dummy, &file, (char*)(dirName+"NumVelocities").c_str());
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void customtwo<T, DIM>::set_default_parameters(){
    double scale = systemInfo->energy_scale;

    NEnergies = 513;
    default_NEnergies = true;

    filename = "custom_two.dat";
    default_filename = true;

    minEnergy = -0.99;
    maxEnergy = 0.99;
    default_energy_limits = true;

    NFermiEnergies = 100;
    minFermiEnergy = -1.0/scale;
    maxFermiEnergy = 1.0/scale;
    default_NFermi = true;
    default_mFermi = true;
    default_MFermi = true;

    deltascat = 0.01/scale;
    scat = 0.01/scale;
    default_scat = true;
    default_deltascat = true;

    temperature = 0.001/scale;
    beta = 1.0/temperature;
    default_temp = true;
}

template <typename T, unsigned DIM>
bool customtwo<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/CustomTwo/";

    get_hdf5(&NumVelocities, &file, (char*)(dirName+"NumVelocities").c_str());

    int numMomentsA, numMomentsB;
    get_hdf5(&numMomentsA, &file, (char*)(dirName+"Vertex0/NumMoments").c_str());
    get_hdf5(&numMomentsB, &file, (char*)(dirName+"Vertex1/NumMoments").c_str());

    if(numMomentsA != numMomentsB){
        file.close();
        return false;
    }
    NumMoments = numMomentsA;

    T temp_from_file;
    try{
        get_hdf5(&temp_from_file, &file, (char*)(dirName+"Temperature").c_str());
        temperature = temp_from_file / systemInfo->energy_scale;
        beta = 1.0/temperature;
        default_temp = false;
    } catch(H5::Exception& e){}

    bool result = false;
    std::string MatrixName = dirName + "Gamma";
    try{
        Gamma = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, NumMoments);
        get_hdf5(Gamma.data(), &file, (char*)MatrixName.c_str());
        result = true;
    } catch(H5::Exception& e) {debug_message("CustomTwo: There is no Gamma matrix.\n");}

    file.close();
    return result;
}

template <typename U, unsigned DIM>
void customtwo<U, DIM>::override_parameters(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    if(variables.CustomTwo_NumEnergies != -1){
        NEnergies = variables.CustomTwo_NumEnergies;
        default_NEnergies = false;
    }

    if(variables.CustomTwo_Temp != -1){
        temperature = variables.CustomTwo_Temp/scale;
        beta = 1.0/temperature;
        default_temp = false;
    }

    if(variables.CustomTwo_Scat != -8888){
        scat = variables.CustomTwo_Scat/scale;
        default_scat = false;
    }

    if(variables.CustomTwo_deltaScat != -8888){
        deltascat = variables.CustomTwo_deltaScat/scale;
        default_deltascat = false;
    } else {
        deltascat = scat;
    }

    if(variables.CustomTwo_FermiMin != -8888){
        minFermiEnergy = (variables.CustomTwo_FermiMin - shift)/scale;
        default_mFermi = false;
    }

    if(variables.CustomTwo_FermiMax != -8888){
        maxFermiEnergy = (variables.CustomTwo_FermiMax - shift)/scale;
        default_MFermi = false;
    }

    if(variables.CustomTwo_NumFermi != -1){
        NFermiEnergies = variables.CustomTwo_NumFermi;
        default_NFermi = false;
    }

    if(variables.CustomTwo_Name != ""){
        filename = variables.CustomTwo_Name;
        default_filename = false;
    }
}

template <typename U, unsigned DIM>
void customtwo<U, DIM>::printCustomTwo(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;
    std::string energy_range = "[" + std::to_string(minEnergy*scale + shift) + ", " + std::to_string(maxEnergy*scale + shift) + "]";

    int p = ((NumVelocities % 4) + 4) % 4;
    std::string component;
    switch(p){
        case 0: component = "-2*Im[Gamma-contraction]"; break;
        case 1: component = "+2*Re[Gamma-contraction]"; break;
        case 2: component = "+2*Im[Gamma-contraction]"; break;
        default: component = "-2*Re[Gamma-contraction]"; break;
    }

    std::cout << "The custom_two() Kubo-Bastin response will be calculated with these parameters: (eV, Kelvin)\n"
        "   Temperature: "             << temperature*scale             << ((default_temp)?         " (default)":"") << "\n"
        "   Broadening: "              << scat*scale                    << ((default_scat)?         " (default)":"") << "\n"
        "   Delta broadening: "        << deltascat*scale               << ((default_deltascat)?    " (default)":"") << "\n"
        "   Max Fermi energy: "        << maxFermiEnergy*scale + shift  << ((default_MFermi)?       " (default)":"") << "\n"
        "   Min Fermi energy: "        << minFermiEnergy*scale + shift  << ((default_mFermi)?       " (default)":"") << "\n"
        "   Number Fermi energies: "   << NFermiEnergies                << ((default_NFermi)?       " (default)":"") << "\n"
        "   Filename: "                << filename                      << ((default_filename)?     " (default)":"") << "\n"
        "   Num integration points: "  << NEnergies                     << ((default_NEnergies)?    " (default)":"") << "\n"
        "   Num Chebyshev moments: "   << NumMoments << "\n"
        "   Number of velocity operators across both vertices: " << NumVelocities <<
        " (NumVelocities mod 4 = " << p << " -> physical component is " << component << ")\n";
}

template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> customtwo<U, DIM>::fill_delta(){
    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR;
    greenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumMoments, NEnergies);

    U factor;
    std::complex<U> complexEnergy;
    for(int i = 0; i < NEnergies; i++){
        complexEnergy = std::complex<U>(energies(i), deltascat);
        for(int m = 0; m < NumMoments; m++){
            factor = -1.0/(1.0 + U(m==0))/M_PI;
            greenR(m, i) = green(m, 1, complexEnergy).imag()*factor;
        }
    }
    return greenR;
}

template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> customtwo<U, DIM>::fill_dgreenR(){
    std::complex<U> complexEnergyP;

    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR;
    dgreenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NEnergies, NumMoments);

    U factor;
    for(int i = 0; i < NEnergies; i++){
        complexEnergyP = std::complex<U>(energies(i), scat);
        for(int m = 0; m < NumMoments; m++){
            factor = 1.0/(1.0 + U(m==0));
            dgreenR(i, m) = dgreen<U>(m, 1, complexEnergyP)*factor;
        }
    }
    return dgreenR;
}

// Unlike conductivity_dc's own triple_product (tools/Src/CondDC/fill.cpp), which
// hardcodes "2*Im[...]" because NumVelocities is always exactly 2 for a built-in
// conductivity, this returns the FULL COMPLEX contraction Z(E) = sum_mn
// delta_m(E) * Gamma_mn * dgreen_n(E) -- the mod-4-dependent component selection
// happens afterward, in calculate(), since it depends on NumVelocities.
template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, 1> customtwo<U, DIM>::triple_product(
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR,
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR){

  Eigen::Matrix<std::complex<U>, -1, -1> GammaEN = Gamma.matrix() * greenR;
  Eigen::Matrix<std::complex<U>, -1, 1> Z = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NEnergies);
  for(int i = 0; i < NEnergies; i++)
      Z(i) = std::complex<U>(dgreenR.row(i) * GammaEN.col(i));
  return Z;
}

template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, 1> customtwo<U, DIM>::calc_response(
  Eigen::Matrix<std::complex<U>, -1, 1> X){

  Eigen::Matrix<std::complex<U>, -1, 1> response =
    Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NFermiEnergies, 1);
  Eigen::Matrix<std::complex<U>, -1, 1> integrand =
    Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NEnergies, 1);

  U fermi;
  for(int i = 0; i < NFermiEnergies; i++){
    fermi = fermiEnergies(i);
    for(int j = 0; j < NEnergies; j++)
      integrand(j) = X(j)*fermi_function(energies(j), fermi, beta);
    response(i) = integrate(energies, integrand);
  }
  return response;
}

template <typename U, unsigned DIM>
void customtwo<U, DIM>::save_to_file(Eigen::Matrix<std::complex<U>, -1, 1> response){
  std::ofstream myfile;
  myfile.open(filename);
  for(int i = 0; i < NFermiEnergies; i++){
    U energy = fermiEnergies(i)*systemInfo->energy_scale + systemInfo->energy_shift;
    myfile << energy << " " << response(i).real() << " " << response(i).imag() << "\n";
  }
  myfile.close();
}

template <typename U, unsigned DIM>
void customtwo<U, DIM>::calculate(){
    if(NEnergies % 2 != 1)
        NEnergies += 1;

    energies = Eigen::Matrix<U, -1, 1>::LinSpaced(NEnergies, minEnergy, maxEnergy);
    fermiEnergies = Eigen::Matrix<U, -1, 1>::LinSpaced(NFermiEnergies, minFermiEnergy, maxFermiEnergy);

    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR = fill_delta();
    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR = fill_dgreenR();

    Eigen::Matrix<std::complex<U>, -1, 1> Z = triple_product(greenR, dgreenR);

    // Derivation (see docs/documentation/examples/custom_vertex_operators.md and
    // rashba_edelstein_graphene_process.py's edelstein() docstring for the full
    // Chebyshev-index bookkeeping): write each raw stored operator as
    // X = i^n_X * X_H with X_H genuinely Hermitian (n_X = count of "missing i"
    // from KITE's raw velocity/commutator token convention). The physical
    // integrand is -2*Im[Tr[A_H delta(E-H) B_H dG+/dE]]; custom_two's stochastic
    // trace gives Gamma_mn = i^p * Tr[T_m(H) A_H T_n(H) B_H], p = NumVelocities.
    // Since i^p has PERIOD 4 (not 2), the correct component of Z depends on
    // p mod 4, not just its parity -- conductivity_dc's own machinery only ever
    // sees p%4==2 (two genuine velocity operators) and hardcodes "2*Im[Z]"
    // accordingly; this generalizes to the other three residues.
    int p = ((NumVelocities % 4) + 4) % 4;
    Eigen::Matrix<U, -1, 1> X_real;
    switch(p){
        case 0: X_real = -2.0*Z.imag(); break;
        case 1: X_real = 2.0*Z.real(); break;
        case 2: X_real = 2.0*Z.imag(); break;
        default: X_real = -2.0*Z.real(); break;
    }
    Eigen::Matrix<std::complex<U>, -1, 1> X = X_real.template cast<std::complex<U>>();

    Eigen::Matrix<std::complex<U>, -1, 1> response = calc_response(X);

    // Normalization validated against the quantized Z2 spin Hall plateau
    // (examples/kane_mele_spin_hall.py, NumVelocities=2 case) -- NOT the same
    // constant as conductivity_dc.cpp's own "den" (that machinery's Gamma comes
    // from a structurally different C++ moment-accumulation path, Gamma2D.cpp,
    // and is independently calibrated against the Haldane Chern-insulator
    // e^2/h plateau; the two are not assumed to share a normalization constant
    // just because both reconstruct a Kubo-Bastin trace).
    U den = 2.0*systemInfo->num_orbitals*systemInfo->spin_degeneracy
            / systemInfo->unit_cell_area / unit_scale;
    U energy_scale_correction = std::pow(systemInfo->energy_scale, NumVelocities - 2);

    // Overall sign: reading Gamma directly into a square Eigen array here (matching
    // conductivity_dc.cpp's own read pattern) amounts to a different Gamma orientation
    // than the Python post-processing scripts use (moments_matrix = f[...][:].T).
    // Empirically verified against two independent, physically-grounded references:
    // for p%4==2 (examples/kane_mele_spin_hall.py), this flips the reconstructed
    // spin Hall plateau's sign relative to the established, documented value
    // (docs/documentation/examples/custom_vertex_operators.md: sigma_sz_xy ~ -2.02,
    // so sign_convention=-1 is needed there); for p%4==1
    // (examples/rashba_edelstein_graphene.py's pure-Rashba case), the orientation
    // difference does NOT flip the sign relative to the from-scratch k-space
    // cross-check (sign_convention=+1 needed there). This is consistent with the
    // even/odd-NumVelocities split derived in rashba_edelstein_graphene_process.py's
    // edelstein() docstring: Im(Gamma) is antisymmetric under m<->n for even
    // NumVelocities (Hermitian Gamma) while Re(Gamma) is antisymmetric for odd
    // NumVelocities (anti-Hermitian Gamma) -- the two cases pick up the Gamma-
    // orientation difference through a different symmetric/antisymmetric channel,
    // so the same fixed "-1" does not apply uniformly across all four residues.
    U sign_convention = (p % 2 == 0) ? -1.0 : 1.0;

    response = response*(sign_convention*den*energy_scale_correction);

    save_to_file(response);
}

template class customtwo<float, 1u>;
template class customtwo<float, 2u>;
template class customtwo<float, 3u>;

template class customtwo<double, 1u>;
template class customtwo<double, 2u>;
template class customtwo<double, 3u>;

template class customtwo<long double, 1u>;
template class customtwo<long double, 2u>;
template class customtwo<long double, 3u>;
