/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/


template <typename T, unsigned DIM>
class singleshot{
    // Single-shot longitudinal DC conductivity.
    //
    // Unlike every other quantity KITE-tools handles, this one needs no
    // Chebyshev reconstruction: KITEx solves the linear system at each energy
    // and stores the converged conductivity itself. All that is left is to put
    // it in a .dat file, which is why this class has no kernel, no broadening
    // and no energy grid of its own -- the energies are the ones that were
    // requested when the calculation was set up.

    public:

        system_info<T, DIM> *systemInfo;    // information about the Hamiltonian
        shell_input variables;              // input from the shell

        bool isRequired;    // was the single-shot conductivity asked for?
        bool isPossible;    // is the data needed to write it in the file?

        int NEnergies;
        int Direction;
        Eigen::Array<double, -1, 1> energies;   // in KPM units
        Eigen::Array<double, -1, 1> broadening; // in KPM units
        Eigen::Array<double, -1, 1> cond;

        std::string filename;
        bool default_filename;

        singleshot(system_info<T, DIM>&, shell_input &);
        singleshot();
        bool is_required();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void printSingleShot();
        void calculate();
};
