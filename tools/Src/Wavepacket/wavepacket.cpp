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
#include "wavepacket.hpp"
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
wavepacket<T, DIM>::wavepacket(system_info<T, DIM>& sysinfo, shell_input & vari){
    H5::Exception::dontPrint();

    systemInfo = &sysinfo;
    variables = vari;

    wavepacket_finished = false;
    isPossible = false;
    isRequired = is_required() && variables.Wavepacket_is_required;
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();
        if(isPossible){
            printWavepacket();
            calculate();
        } else {
            std::cout << "ERROR. gaussian_wave_packet() post-processing was requested but "
                "the data needed for it was not found in the input .h5 file. Note KITEx "
                "only builds the wave-packet driver when it was compiled with a gcc new "
                "enough to provide the Bessel functions. Exiting.";
            exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool wavepacket<T, DIM>::is_required(){
    std::string name = systemInfo->filename;
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    bool result = false;
    try{
        get_hdf5(&NumDisorder, &file, (char*)"/Calculation/gaussian_wave_packet/NumDisorder");
        result = true;
    } catch(H5::Exception& e){}

    file.close();
    return result;
}

template <typename T, unsigned DIM>
void wavepacket<T, DIM>::set_default_parameters(){
    filename = "wavepacket.dat";
    default_filename = true;
}

template <typename T, unsigned DIM>
void wavepacket<T, DIM>::override_parameters(){
    if(variables.Wavepacket_Name != ""){
        filename         = variables.Wavepacket_Name;
        default_filename = false;
    }
}

template <typename T, unsigned DIM>
void wavepacket<T, DIM>::printWavepacket(){
    const double scale = systemInfo->energy_scale;

    std::cout << "The gaussian_wave_packet() time series will be written with these parameters:\n"
        "   Time points: "         << NumPoints   << "\n"
        "   Timestep: "            << timestep << " (rescaled), "
                                   << timestep/scale << " hbar/eV\n"
        "   Chebyshev moments per step: " << NumMoments  << "\n"
        "   Disorder realisations: "      << NumDisorder << "\n"
        "   Packet width: "        << width << "\n"
        "   Tracked operators: "   << operator_labels.size() << "\n"
        "   Filename: "            << filename << ((default_filename)? " (default)":"") << "\n";
}

template <typename T, unsigned DIM>
bool wavepacket<T, DIM>::fetch_parameters(){
    std::string name = systemInfo->filename;
    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    const std::string dirName = "/Calculation/gaussian_wave_packet/";

    get_hdf5(&NumDisorder, &file, (char*)(dirName+"NumDisorder").c_str());
    get_hdf5(&NumMoments,  &file, (char*)(dirName+"NumMoments").c_str());
    get_hdf5(&NumPoints,   &file, (char*)(dirName+"NumPoints").c_str());
    get_hdf5(&width,       &file, (char*)(dirName+"width").c_str());

    // timestep is stored single-precision whatever the run's precision is.
    float stored_timestep;
    get_hdf5(&stored_timestep, &file, (char*)(dirName+"timestep").c_str());
    timestep = stored_timestep;

    // Every series is written as an Eigen (NumPoints, 1) complex array, one per
    // axis and one per operator label. Id is not read: see wavepacket.hpp.
    const std::string axes = "xyz";
    bool result = false;
    try{
        for(unsigned d = 0; d < DIM; d++){
            Eigen::Array<std::complex<T>,-1,-1> m =
                Eigen::Array<std::complex<T>,-1,-1>::Zero(NumPoints, 1);
            Eigen::Array<std::complex<T>,-1,-1> v = m;
            get_hdf5(m.data(), &file, (char*)(dirName+"mean_value"+axes[d]).c_str());
            get_hdf5(v.data(), &file, (char*)(dirName+"Var"+axes[d]).c_str());
            mean_value.push_back(m);
            variance.push_back(v);
        }
        result = true;
    } catch(H5::Exception& e) {debug_message("wave packet: there is no time series.\n");}

    try{
        read_label_list(operator_labels, file, dirName + "Operators");
        for(const auto &label : operator_labels){
            Eigen::Array<std::complex<T>,-1,-1> s =
                Eigen::Array<std::complex<T>,-1,-1>::Zero(NumPoints, 1);
            get_hdf5(s.data(), &file, (char*)(dirName+label).c_str());
            operator_series.push_back(s);
        }
    } catch(H5::Exception& e) {
        debug_message("wave packet: no operators were tracked.\n");
        operator_labels.clear();
        operator_series.clear();
    }

    file.close();
    return result;
}

template <typename U, unsigned DIM>
void wavepacket<U, DIM>::calculate(){
    const double scale = systemInfo->energy_scale;
    const std::string axes = "xyz";

    std::ofstream myfile;
    myfile.open(filename);
    myfile << "# Wave-packet time evolution, from calculation.gaussian_wave_packet().\n";
    myfile << "# Averaged over " << NumDisorder << " disorder realisation"
           << ((NumDisorder == 1)? "" : "s") << ", "
           << NumMoments << " Chebyshev moments per step.\n";
    myfile << "# KITEx propagates the rescaled Hamiltonian, so the stored timestep is the\n"
              "# argument of the Bessel expansion of exp(-i*Htilde*ttilde). Column 1 is the\n"
              "# physical time t = t_rescaled/EnergyScale, in hbar/eV with EnergyScale = "
           << scale << " eV;\n"
              "# column 2 is t_rescaled exactly as stored.\n";
    myfile << "# Positions and variances are Cartesian, carrying the length unit and its\n"
              "# square from the lattice vectors supplied in the model. The operator columns\n"
              "# carry whatever units the registered vertex was given.\n";

    myfile << "# col 1: t   col 2: t_rescaled";
    int col = 2;
    for(unsigned d = 0; d < DIM; d++)
        myfile << "   col " << ++col << ": <" << axes[d] << ">";
    for(unsigned d = 0; d < DIM; d++)
        myfile << "   col " << ++col << ": Var" << axes[d];
    for(const auto &label : operator_labels)
        myfile << "   col " << ++col << ": <" << label << ">";
    myfile << "\n";

    for(int t = 0; t < NumPoints; t++){
        const double trescaled = t*timestep;
        myfile << trescaled/scale << " " << trescaled;
        for(unsigned d = 0; d < DIM; d++)
            myfile << " " << std::real(mean_value[d](t, 0));
        for(unsigned d = 0; d < DIM; d++)
            myfile << " " << std::real(variance[d](t, 0));
        for(std::size_t k = 0; k < operator_series.size(); k++)
            myfile << " " << std::real(operator_series[k](t, 0));
        myfile << "\n";
    }
    myfile.close();
    wavepacket_finished = true;
}

template class wavepacket<float, 1u>;
template class wavepacket<float, 2u>;
template class wavepacket<float, 3u>;

template class wavepacket<double, 1u>;
template class wavepacket<double, 2u>;
template class wavepacket<double, 3u>;

template class wavepacket<long double, 1u>;
template class wavepacket<long double, 2u>;
template class wavepacket<long double, 3u>;
