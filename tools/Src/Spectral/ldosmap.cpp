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
#include "ldosmap.hpp"
#include "../Tools/functions.hpp"

#include "../macros.hpp"

// See ldos.cpp for why the variable-length string list is read by hand.
static void read_label_list(std::vector<std::string> &labels, H5::H5File &file,
                            const std::string &path){
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
ldosmap<T, DIM>::ldosmap(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    ldosmap_finished = false;
    isPossible = false;
    isRequired = is_required() && variables.LDOSmap_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        if(isPossible){
            printLDOSmap();
            calculate();
        } else {
            std::cout << "ERROR. ldos_map() post-processing was requested but the data "
                "needed for it (/Calculation/ldos_map/Map) was not found in the input .h5 "
                "file. Note KITEx only builds the map drivers for a complex Hamiltonian: "
                "with is_complex=False it reports success and writes no map. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool ldosmap<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    bool result = false;
    try{
        get_hdf5(&NumVectors, &file, (char*)"/Calculation/ldos_map/NumVectors");
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void ldosmap<T, DIM>::set_default_parameters(){
    filename = "ldos_map.dat";
    default_filename = true;
}

template <typename T, unsigned DIM>
void ldosmap<T, DIM>::override_parameters(){
    if(variables.LDOSmap_Name != ""){
        filename         = variables.LDOSmap_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void ldosmap<T, DIM>::printLDOSmap(){
    double scale = systemInfo->energy_scale;
    double shift = systemInfo->energy_shift;

    std::cout << "The ldos_map() real-space map will be written with these parameters: (eV)\n"
        "   Energy: "              << Energy*scale + shift << "\n"
        "   Broadening (sigma): "  << Sigma*scale << "\n"
        "   Number of vectors: "   << NumVectors  << "\n"
        "   Unit cells: "          << Nt << ", orbitals: " << Orb << "\n"
        "   Coefficients: "        << ((Coef_ID)? "window (dimensionless)"
                                                : "gaussian (1/eV)") << "\n"
        "   Filename: "            << filename << ((default_filename)? " (default)":"") << "\n";
}

template <typename T, unsigned DIM>
bool ldosmap<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/ldos_map/";

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
    } catch(H5::Exception& e) {debug_message("ldos_map: There is no Map.\n");}

    // Optional operator-weighted maps, one value per unit cell (summed over
    // orbitals), written only when calculation.ldos_map(operators=[...]) was used.
    try{
        read_label_list(operator_labels, file, dirName + "Operators");
        for(const auto &label : operator_labels){
            Eigen::Array<T,-1,-1> m = Eigen::Array<T,-1,-1>::Zero(Nt, 2);
            get_hdf5(m.data(), &file, (char*)(dirName+"Map_Operators/"+label).c_str());
            Map_Operators.push_back(m);
        }
    } catch(H5::Exception& e) {
        debug_message("ldos_map: no operator-weighted maps requested.\n");
        operator_labels.clear();
        Map_Operators.clear();
    }

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void ldosmap<T, DIM>::write_one(const Eigen::Array<T, -1, -1> &data, long entries,
                                const std::string &out, const std::string &what){
    const double scale = systemInfo->energy_scale;
    const double shift = systemInfo->energy_shift;

    std::ofstream myfile;
    myfile.open(out);
    myfile << "# " << what << ", from calculation.ldos_map(), at E = "
           << Energy*scale + shift << " eV, sigma = " << Sigma*scale << " eV.\n";
    if(Coef_ID)
        myfile << "# Coefficients: window. Dimensionless -- the spectral weight in a\n"
                  "# window of width sigma about E, not a density.\n";
    else
        myfile << "# Coefficients: gaussian. Units 1/eV -- a density of states per site,\n"
                  "# Gaussian-broadened with standard deviation sigma in real energy units.\n";
    myfile << "# Averaged over " << NumVectors << " random vectors; the last column is the\n"
              "# standard error of that mean.\n";
    if(DIM == 2)
        myfile << "# col 1: x   col 2: y   col 3: orbital   col 4: value   col 5: error\n";
    else
        myfile << "# col 1: x   col 2: y   col 3: z   col 4: orbital   col 5: value   col 6: error\n";

    const long Lx = systemInfo->size[0];
    const long Ly = (DIM > 1) ? systemInfo->size[1] : 1;
    // index = x + y*Lx [+ z*Lx*Ly] + orbital*Nt; x fastest, orbital slowest.
    for(long i = 0; i < entries; i++){
        const long orbital = (entries == Sizet) ? i / Nt : 0;
        const long cell = i - orbital*Nt;
        const long x = cell % Lx;
        const long y = (cell / Lx) % Ly;
        myfile << x << " " << y << " ";
        if(DIM == 3) myfile << cell / (Lx*Ly) << " ";
        if(entries == Sizet) myfile << orbital << " ";
        else                 myfile << "all ";
        myfile << std::real(data(i, 0)) << " " << std::real(data(i, 1)) << "\n";
    }
    myfile.close();
}

template <typename U, unsigned DIM>
void ldosmap<U, DIM>::calculate(){
    write_one(Map, Sizet, filename, "Local density of states");

    // One file per operator label, named after it, so several labels from one
    // run do not collide.
    for(std::size_t k = 0; k < operator_labels.size(); k++){
        std::string out = filename;
        const std::size_t dot = out.rfind('.');
        const std::string stem = (dot == std::string::npos) ? out : out.substr(0, dot);
        const std::string ext  = (dot == std::string::npos) ? "" : out.substr(dot);
        write_one(Map_Operators[k], Nt, stem + "_" + operator_labels[k] + ext,
                  "Operator-weighted local density of states (" + operator_labels[k] +
                  "), summed over orbitals");
    }
    ldosmap_finished = true;
}

template class ldosmap<float, 1u>;
template class ldosmap<float, 2u>;
template class ldosmap<float, 3u>;

template class ldosmap<double, 1u>;
template class ldosmap<double, 2u>;
template class ldosmap<double, 3u>;

template class ldosmap<long double, 1u>;
template class ldosmap<long double, 2u>;
template class ldosmap<long double, 3u>;
