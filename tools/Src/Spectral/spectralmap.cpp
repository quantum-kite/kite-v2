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
#include "spectralmap.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
spectralmap<T, DIM>::spectralmap(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    spectralmap_finished = false;
    isPossible = false;
    isRequired = is_required() && variables.SpectralMap_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        if(isPossible){
            printSpectralMap();
            calculate();
        } else {
            std::cout << "ERROR. spectral_map() post-processing was requested but the data "
                "needed for it (/Calculation/spectral_map/Map) was not found in the input .h5 "
                "file. Note KITEx only builds the map drivers for a complex Hamiltonian: "
                "with is_complex=False it reports success and writes no map. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool spectralmap<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    bool result = false;
    try{
        get_hdf5(&NumVectors, &file, (char*)"/Calculation/spectral_map/NumVectors");
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void spectralmap<T, DIM>::set_default_parameters(){
    filename = "spectral_map.dat";
    default_filename = true;
}

template <typename T, unsigned DIM>
void spectralmap<T, DIM>::override_parameters(){
    if(variables.SpectralMap_Name != ""){
        filename         = variables.SpectralMap_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void spectralmap<T, DIM>::printSpectralMap(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    std::cout << "The spectral_map() momentum-space map will be written with these parameters: (eV)\n"
        "   Energy: "              << Energy*scale + shift << "\n"
        "   Broadening (sigma): "  << Sigma*scale << "\n"
        "   Number of vectors: "   << NumVectors  << "\n"
        "   k-points: "            << Nt << ", orbitals: " << Orb << "\n"
        "   Coefficients: "        << ((Coef_ID)? "window (dimensionless)"
                                                : "gaussian (1/eV)") << "\n"
        "   Filename: "            << filename << ((default_filename)? " (default)":"") << "\n";
}

template <typename T, unsigned DIM>
bool spectralmap<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/spectral_map/";

    get_hdf5(&NumVectors, &file, (char*)(dirName+"NumVectors").c_str());
    get_hdf5(&Coef_ID,    &file, (char*)(dirName+"Coef_ID").c_str());
    get_hdf5(&Energy,     &file, (char*)(dirName+"Energy").c_str());
    get_hdf5(&Sigma,      &file, (char*)(dirName+"Sigma").c_str());

    Orb = systemInfo->num_orbitals;
    Nt = 1;
    for(unsigned d = 0; d < DIM; d++)
        Nt *= systemInfo->size[d];
    Sizet = Nt * Orb;

    // Stored as (2, N) row-major, i.e. an Eigen (N, 2) column-major array:
    // column 0 the value, column 1 the standard error over NumVectors draws.
    bool result = false;
    try{
        Map = Eigen::Array<T,-1,-1>::Zero(Sizet, 2);
        get_hdf5(Map.data(), &file, (char*)(dirName+"Map").c_str());
        result = true;
    } catch(H5::Exception& e) {debug_message("spectral_map: There is no Map.\n");}

    file.close();
    return result;
}

template <typename U, unsigned DIM>
void spectralmap<U, DIM>::calculate(){
    const double scale = systemInfo->energy_scale;
    const double shift = systemInfo->energy_shift;

    // Same conversion arpes.cpp uses, so the two readers quote k in the same
    // frame: a k given in fractional reciprocal-lattice coordinates is
    // k_cartesian = 2*pi * vectors^{-T} * k_fractional.
    const Eigen::Matrix<double, -1, -1> recip =
        2.0 * M_PI * systemInfo->vectors.matrix().inverse().transpose();

    std::ofstream myfile;
    myfile.open(filename);
    myfile << "# Momentum-resolved spectral function, from calculation.spectral_map(), at E = "
           << Energy*scale + shift << " eV, sigma = " << Sigma*scale << " eV.\n";
    if(Coef_ID)
        myfile << "# Coefficients: window. Dimensionless -- approximately the spectral weight\n"
                  "# in a window of width sigma centred on E. The estimator is W^2, so the\n"
                  "# Jackson smoothing at the window edges makes it a slight underestimate.\n";
    else
        myfile << "# Coefficients: gaussian. Units 1/eV -- the momentum-resolved spectral\n"
                  "# function A(k, orbital; E), diagonal in orbital, Gaussian-broadened with\n"
                  "# standard deviation sigma in real energy units. Integrating over E gives\n"
                  "# 1 for each (k, orbital).\n";
    myfile << "# Averaged over " << NumVectors << " random vectors; the last column is the\n"
              "# standard error of that mean, over random vectors (and twists) only. The\n"
              "# disorder is drawn once, so it carries no realization-to-realization spread.\n";
    myfile << "# k is given twice: first in fractional reciprocal-lattice coordinates f, in\n"
              "# [-1/2, 1/2), with k.a_d = 2*pi*f_d -- the 2*pi is in the reciprocal vectors,\n"
              "# not in f. That range is one full Brillouin zone as a parallelepiped in the\n"
              "# b_d basis, not the Wigner-Seitz zone, so on a honeycomb lattice the sampled\n"
              "# region is a rhombus rather than a hexagon. Then the same point in Cartesian\n"
              "# components, carrying the inverse of whatever length unit the lattice vectors\n"
              "# were given in.\n";
    myfile << "# The k grid is fixed by the simulation cell, so these labels are only\n"
              "# meaningful under periodic or fixed-twist boundaries. With random boundary\n"
              "# twists each vector is drawn at a different twist and the average smears\n"
              "# neighbouring k cells together.\n";
    if(DIM == 1)
        myfile << "# col 1: k1   col 2: kx   col 3: orbital   col 4: value   col 5: error\n";
    else if(DIM == 2)
        myfile << "# col 1: k1   col 2: k2   col 3: kx   col 4: ky   "
                  "col 5: orbital   col 6: value   col 7: error\n";
    else
        myfile << "# col 1: k1   col 2: k2   col 3: k3   col 4: kx   col 5: ky   col 6: kz   "
                  "col 7: orbital   col 8: value   col 9: error\n";

    const long Lx = systemInfo->size[0];
    const long Ly = (DIM > 1) ? systemInfo->size[1] : 1;

    Eigen::Matrix<double, -1, 1> frac(DIM);
    for(long i = 0; i < Sizet; i++){
        const long orbital = i / Nt;
        const long cell = i - orbital*Nt;
        const long s[3] = {cell % Lx, (cell / Lx) % Ly, cell / (Lx*Ly)};

        // Undo the (s + L/2) % L re-indexing store_spectral applied: index
        // L_d/2 is k_d = 0, and the fractional coordinate runs over [-1/2, 1/2).
        for(unsigned d = 0; d < DIM; d++){
            // Integer L/2, to match store_spectral's own integer
            // (coord + N/2) % N. KITEx refuses a length that is not a multiple
            // of TILE, so L is even and the two agree; this just keeps the
            // reader exact rather than exact-by-coincidence.
            const long L = systemInfo->size[d];
            frac(d) = double(s[d] - L/2)/L;
            myfile << frac(d) << " ";
        }
        const Eigen::Matrix<double, -1, 1> kcart = recip * frac;
        for(unsigned d = 0; d < DIM; d++)
            myfile << kcart(d) << " ";

        myfile << orbital << " "
               << std::real(Map(i, 0)) << " " << std::real(Map(i, 1)) << "\n";
    }
    myfile.close();
    spectralmap_finished = true;
}

template class spectralmap<float, 1u>;
template class spectralmap<float, 2u>;
template class spectralmap<float, 3u>;

template class spectralmap<double, 1u>;
template class spectralmap<double, 2u>;
template class spectralmap<double, 3u>;

template class spectralmap<long double, 1u>;
template class spectralmap<long double, 2u>;
template class spectralmap<long double, 3u>;
