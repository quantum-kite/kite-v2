/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Reconstructs the energy-resolved spectral density of a custom_one() rank-one
// vertex A, Tr[A * delta(E-H)], from its raw Chebyshev moments -- the exact same
// Jackson-kernel Chebyshev-to-density machinery as dos<T,DIM> (Spectral/dos.hpp),
// just applied to Gamma (the vertex-weighted moments) instead of MU (the trivial,
// identity-operator DOS moments). See customone.cpp for why NumVelocities (auto-
// detected from the vertex string on the Python side, src/kite/__init__.py) matters:
// KITE's raw "v" token is missing a factor of i, so Tr[T_n(H)*A] comes out real for
// an even total count of velocity operators in A, purely imaginary for an odd count
// -- both are the genuine physical signal, just in a different component.

template <typename T, unsigned DIM>
class customone{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        double Emax, Emin;
        bool default_Emax, default_Emin;

        int NumMoments;
        int MaxMoments;
        bool default_NumMoments;
        int NumVelocities;

        std::string kernel;
        double kernel_parameter;
        bool default_kernel;
        bool default_kernel_parameter;

        std::string filename;
        bool default_filename;

        bool default_NEnergies;
        int NEnergies;
        Eigen::Matrix<T, -1, 1> energies;

        Eigen::Array<std::complex<T>, -1, -1> Gamma;
        Eigen::Array<std::complex<T>, -1, -1> GammaE;
        bool customone_finished;

        customone(system_info<T, DIM>&, shell_input &);
        customone();
        bool is_required();
        void printCustomOne();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();

};
