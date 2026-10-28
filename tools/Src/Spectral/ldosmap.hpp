/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Decodes the real-space maps written by calculation.ldos_map() into a .dat
// with site coordinates, so a user does not have to reimplement the index
// arithmetic in Python -- three shipped processors do it by hand today.
//
// There is no reconstruction to perform. Unlike every other unit here, KITEx
// finishes the energy resolution itself: calc_ldos sums the Chebyshev series
// into a single vector against a Gaussian (or window) coefficient set, and
// stores the resulting density per site. This unit only decodes the flat site
// index and converts the energy back to eV.
//
// Layout, verified against a run rather than read off: a vacancy returns an
// LDoS of exactly zero, and on a 16x24 lattice a vacancy at (x=5, y=3) on the
// first orbital zeroes column 53 = 5 + 3*16 + 0*384 and no other. So the flat
// index is
//
//     index = x + y*Lx [+ z*Lx*Ly] + orbital*Nt,     Nt = Lx*Ly[*Lz]
//
// x fastest, orbital slowest, and NOT fftshifted. The operator variant carries
// one value per site summed over orbitals, so it has Nt entries, not Sizet.

template <typename T, unsigned DIM>
class ldosmap{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        int NumVectors;
        int Coef_ID;
        double Energy;
        double Sigma;

        long Nt;          // unit cells
        long Sizet;       // unit cells * orbitals
        int Orb;

        std::string filename;
        bool default_filename;

        Eigen::Array<T, -1, -1> Map;                      // (Sizet, 2): value, error
        std::vector<std::string> operator_labels;
        std::vector<Eigen::Array<T, -1, -1>> Map_Operators;

        bool ldosmap_finished;

        ldosmap(system_info<T, DIM>&, shell_input &);
        ldosmap();
        bool is_required();
        void printLDOSmap();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();
        void write_one(const Eigen::Array<T, -1, -1> &, long, const std::string &,
                       const std::string &);
};
