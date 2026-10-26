/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Reconstructs the local, operator-weighted spectral density of a
// custom_one_local() vertex A at each probed orbital, from the raw Chebyshev
// moments Gamma_n^(j) = <j|T_n(Htilde)*A|j> stored at
// /Calculation/CustomOneLocal/Gamma.
//
// Same Jackson-kernel machinery as customone<T,DIM>, with two differences that
// follow from the seeding being a single basis ket rather than a random vector:
//
//   * no 1/N -- there is no stochastic trace to normalize;
//   * one output column per probed orbital, instead of one curve.
//
// store_custom_one_local (Src/Simulation/Custom/SimulationRankOneLocal.cpp)
// stores 0.5*(Gamma + conj(Gamma)), so what reaches the file is already the
// symmetrized combination <j|{T_n,A}|j>/2, which is real. Unlike customone, no
// velocity-parity rotation is applied or possible: Python writes no
// NumVelocities for the Local variants, and the single-site matrix element has
// no trace cyclicity to make an anti-Hermitian A's moments purely imaginary.
// For a vertex that IS diagonal in the orbital basis -- the intended on-site
// use -- A|j> = A_jj|j>, so Gamma_n = A_jj*mu_n exactly and the reconstruction
// is the orbital-resolved LDoS scaled by A_jj.

template <typename T, unsigned DIM>
class customonelocal{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        double Emax, Emin;
        bool default_Emax, default_Emin;

        int NumMoments;
        int MaxMoments;
        long NumPositions;

        std::string filename;
        bool default_filename;

        bool default_NEnergies;
        int NEnergies;
        Eigen::Matrix<T, -1, 1> energies;

        Eigen::Matrix<unsigned long, -1, -1> Orbitals;
        Eigen::Matrix<unsigned long, -1, -1> Positions;
        Eigen::Array<std::complex<T>, -1, -1> Gamma;
        Eigen::Array<std::complex<T>, -1, -1> GammaE;
        bool customonelocal_finished;

        customonelocal(system_info<T, DIM>&, shell_input &);
        customonelocal();
        bool is_required();
        void printCustomOneLocal();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();

};
