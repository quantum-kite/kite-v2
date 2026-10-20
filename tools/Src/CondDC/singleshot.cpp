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
#include "H5Cpp.h"
#include "../Tools/ComplexTraits.hpp"
#include "../Tools/myHDF5.hpp"

#include "../Tools/parse_input.hpp"
#include "../Tools/systemInfo.hpp"
#include "singleshot.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
singleshot<T, DIM>::singleshot(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    isPossible = false;
    isRequired = is_required() && variables.SingleShot_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        if(isPossible){
            printSingleShot();
            calculate();
        } else {
            std::cout << "ERROR. The single-shot DC conductivity was requested "
                "but the data needed to write it was not found in the input .h5 "
                "file. Make sure KITEx has processed the file first. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool singleshot<T, DIM>::is_required(){
    // Was a single-shot calculation set up in this file? The group exists only
    // if calculation.singleshot_conductivity_dc() was called.

    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    bool result = false;
    try{
        get_hdf5(&Direction, &file,
                 (char*)"/Calculation/singleshot_conductivity_dc/Direction");
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void singleshot<T, DIM>::set_default_parameters(){
    filename = "singleshot_condDC.dat";
    default_filename = true;
}

template <typename T, unsigned DIM>
bool singleshot<T, DIM>::fetch_parameters(){
    // Reads the result KITEx already computed. Each row of `SingleShot` holds
    // energy and broadening in eV, the preserve_disorder flag, and the
    // conductivity. Only the last of those is taken from it -- the energies are
    // read from `Energy` instead, in KPM units, so that the shift back to eV
    // here is the same one every other KITE-tools quantity applies.

    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    std::string dirName = "/Calculation/singleshot_conductivity_dc/";

    bool result = false;
    try{
        const H5::DataSet dataset = file.openDataSet(dirName + "Energy");
        const H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dim[2];
        dataspace.getSimpleExtentDims(dim, NULL);
        NEnergies = int(dim[0]*dim[1]);

        energies   = Eigen::Array<double, -1, 1>::Zero(NEnergies);
        broadening = Eigen::Array<double, -1, 1>::Zero(NEnergies);
        get_hdf5<double>(energies.data(), &file, (char*)(dirName+"Energy").c_str());
        get_hdf5<double>(broadening.data(), &file, (char*)(dirName+"Gamma").c_str());

        // Stored as (4, NEnergies) by SimulationSingleShot.cpp, so the four
        // entries belonging to one energy are contiguous.
        Eigen::Array<double, -1, -1> table =
            Eigen::Array<double, -1, -1>::Zero(4, NEnergies);
        get_hdf5<double>(table.data(), &file, (char*)(dirName+"SingleShot").c_str());
        cond = table.row(3).transpose();

        result = true;
    } catch(H5::Exception& e){
        debug_message("SingleShot: no stored conductivity.\n");
    }

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void singleshot<T, DIM>::override_parameters(){
    if(variables.SingleShot_Name != ""){
        filename = variables.SingleShot_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void singleshot<T, DIM>::printSingleShot(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;
    std::string dir = "xx";
    if(Direction == 1) dir = "yy";
    if(Direction == 2) dir = "zz";

    verbose_message("The program will calculate the single-shot DC conductivity with the following parameters:\n");
    verbose_message("   direction: "); verbose_message(dir); verbose_message("\n");
    verbose_message("   number of energies: "); verbose_message(NEnergies); verbose_message("\n");
    verbose_message("   energy range: ["); verbose_message(energies.minCoeff()*scale + shift);
    verbose_message(", "); verbose_message(energies.maxCoeff()*scale + shift); verbose_message("] eV\n");
    verbose_message("   broadening: "); verbose_message(broadening.minCoeff()*scale);
    verbose_message(" to "); verbose_message(broadening.maxCoeff()*scale); verbose_message(" eV\n");
    verbose_message("   file name: "); verbose_message(filename); verbose_message("\n");
    verbose_message("The single-shot method solves for each energy directly, so there is nothing to reconstruct here.\n");
}

template <typename T, unsigned DIM>
void singleshot<T, DIM>::calculate(){
    // Nothing is computed: the conductivity KITEx stored is the final number.
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    std::ofstream myfile;
    myfile.open(filename);
    for(int i = 0; i < NEnergies; i++){
        myfile << energies(i)*scale + shift << " " << cond(i) << "\n";
    }
    myfile.close();
}


// Instantiations
template class singleshot<float, 1u>;
template class singleshot<float, 2u>;
template class singleshot<float, 3u>;

template class singleshot<double, 1u>;
template class singleshot<double, 2u>;
template class singleshot<double, 3u>;

template class singleshot<long double, 1u>;
template class singleshot<long double, 2u>;
template class singleshot<long double, 3u>;
