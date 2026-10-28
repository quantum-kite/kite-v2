/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Decodes the momentum-space maps written by calculation.spectral_map() into a
// .dat with the k-point of every row spelled out, so a user does not have to
// reimplement the index arithmetic -- the shipped Weyl processor does it by
// hand today.
//
// As with ldos_map there is no reconstruction to perform: KITEx finishes the
// energy resolution itself and stores A(k) at a single energy. This unit only
// decodes the flat index and converts the energy back to eV.
//
// Layout, verified against a run rather than read off. On a 16x32 lattice of
// nearly-decoupled chains (t_x = 1, t_y = 0.05, so E(k) depends essentially
// only on k_x) the weight at E = -2cos(pi/4) sits at exactly two columns,
// s_x = 6 and s_x = 10, for every s_y; with the hoppings exchanged it sits at
// s_y = 8 and s_y = 24 for every s_x. So the flat index is
//
//     index = s_x + s_y*Lx [+ s_z*Lx*Ly] + orbital*Nt,   Nt = Lx*Ly[*Lz]
//
// x fastest and orbital slowest, the same as ldos_map -- but unlike ldos_map
// every axis IS fftshifted by store_spectral, so s_d = L_d/2 is k_d = 0 and
//
//     k_d = 2*pi*(s_d - L_d/2)/L_d.
//
// The energy was chosen to make that second statement testable: k = +-pi/2 is
// a fixed point of the (s + L/2) % L re-indexing, so a contour sitting there
// cannot tell a shifted layout from an unshifted one. pi/4 can, and does --
// an unshifted layout would have put the weight at s_x = 2 and s_x = 14.

template <typename T, unsigned DIM>
class spectralmap{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        int NumVectors;
        int Coef_ID;
        double Energy;
        double Sigma;

        long Nt;          // unit cells, i.e. k-points
        long Sizet;       // k-points * orbitals
        int Orb;

        std::string filename;
        bool default_filename;

        Eigen::Array<T, -1, -1> Map;                      // (Sizet, 2): value, error

        bool spectralmap_finished;

        spectralmap(system_info<T, DIM>&, shell_input &);
        spectralmap();
        bool is_required();
        void printSpectralMap();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();
};
