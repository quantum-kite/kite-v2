/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Reconstructs the Kubo-Bastin energy integral of a custom_two() rank-two
// vertex pair [A, B], Gamma_mn = Tr[T_m(H) A T_n(H) B], as a function of
// Fermi energy -- the general-vertex-pair analogue of conductivity_dc's own
// Gamma2D/CondDC machinery (tools/Src/CondDC/), which only ever needs to
// handle the single case of two genuine velocity operators (NumVelocities=2).
//
// custom_two() allows ANY vertex pair (bare density operators, position
// operators, arbitrary orbital-coupling matrices, mixed with velocity
// tokens), so the combined "missing factor of i" bookkeeping (see
// custom_one's own docstring/customone.hpp) needs the FULL mod-4 treatment,
// not conductivity_dc's hardcoded mod-2-equivalent "always take 2*Im[...]"
// shortcut, which is only correct when NumVelocities%4==2 (true for
// conductivity_dc and for symmetrized-current vertex pairs like spin Hall's
// (1/2){v_x,s_z}, but not for e.g. a bare-density-operator vertex like the
// Rashba-Edelstein effect's A=s_x, NumVelocities=1).
//
// See customtwo.cpp for the full derivation.

template <typename T, unsigned DIM>
class customtwo{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        int NumMoments;
        int NumVelocities;

        double temperature;
        T beta;
        bool default_temp;

        T scat, deltascat;
        bool default_scat, default_deltascat;

        int NEnergies;
        bool default_NEnergies;
        T minEnergy, maxEnergy;
        bool default_energy_limits;
        Eigen::Matrix<T, -1, 1> energies;

        int NFermiEnergies;
        double minFermiEnergy, maxFermiEnergy;
        Eigen::Matrix<T, -1, 1> fermiEnergies;
        bool default_NFermi, default_mFermi, default_MFermi;

        std::string filename;
        bool default_filename;

        // As stored by store_custom_two (Src/Simulation/Custom/SimulationRankTwo.cpp):
        // shape (NumMoments, NumMoments), same read convention as conductivity_dc's
        // own Gamma (a genuinely square matrix requires equal moment counts on both
        // vertices -- store_custom_two only Hermitizes, and this reconstruction only
        // applies, when that holds; unequal per-vertex moment counts are rejected).
        Eigen::Array<std::complex<T>, -1, -1> Gamma;

        customtwo(system_info<T, DIM>&, shell_input &);
        bool is_required();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void printCustomTwo();
        void calculate();

        Eigen::Matrix<std::complex<T>, -1, -1, Eigen::ColMajor> fill_delta();
        Eigen::Matrix<std::complex<T>, -1, -1, Eigen::RowMajor> fill_dgreenR();
        Eigen::Matrix<std::complex<T>, -1, 1> triple_product(
          Eigen::Matrix<std::complex<T>, -1, -1, Eigen::ColMajor>,
          Eigen::Matrix<std::complex<T>, -1, -1, Eigen::RowMajor>);
        Eigen::Matrix<std::complex<T>, -1, 1> calc_response(
          Eigen::Matrix<std::complex<T>, -1, 1>);
        void save_to_file(Eigen::Matrix<std::complex<T>, -1, 1>);
};
