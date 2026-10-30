/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2026, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

// Decodes the time series written by calculation.gaussian_wave_packet() into a
// .dat with one row per timestep: the packet's mean position and variance along
// each axis, and the expectation value of every operator the user registered.
//
// There is nothing to reconstruct. KITEx runs the Bessel-series propagator and
// stores the finished observables, one (1, NumPoints) dataset each.
//
// Two things this unit must get right, neither of them arithmetic:
//
// 1. It does NOT emit `Id`. `Id` is the amplitude at a probing point whose
//    location depends on how the lattice was split across threads, so the same
//    physics gives a different `Id` at a different thread count. It is hidden,
//    not explained (AC-WAVE-06), and a reader that simply dumped every dataset
//    in the group would put it back in front of users.
//
// 2. The time axis is stored RESCALED. `timestep` is the argument of the Bessel
//    expansion of exp(-i*Htilde*ttilde), so a step of dt in physical time is
//    dt*EnergyScale there. Both axes are written: `t` in hbar/eV, which is what
//    `python/kite/__init__.py`'s own docstring tells a user to recover by
//    dividing by EnergyScale, and `t_rescaled` exactly as stored.
//
// Layout, verified against a run rather than read off: the spin-precession
// example seeds a spin-1/2 in a uniform transverse field, where <Sz(t)> is
// exactly (1/2)cos(2Bt) and <Sy(t)> is -(1/2)sin(2Bt). The registered `l2` and
// `l1` columns match those to correlation 1.000000 -- on the PHYSICAL axis
// t = n*timestep/EnergyScale. On the rescaled axis the period is wrong by the
// factor EnergyScale, which is 5 in that example.

template <typename T, unsigned DIM>
class wavepacket{
	public:

        system_info<T, DIM> *systemInfo;
        shell_input variables;

        bool isRequired;
        bool isPossible;

        int NumPoints;
        int NumMoments;
        int NumDisorder;
        double timestep;        // rescaled
        double width;

        std::string filename;
        bool default_filename;

        // (NumPoints, 1) each. Id is deliberately not among them.
        std::vector<Eigen::Array<std::complex<T>, -1, -1>> mean_value;
        std::vector<Eigen::Array<std::complex<T>, -1, -1>> variance;
        std::vector<std::string> operator_labels;
        std::vector<Eigen::Array<std::complex<T>, -1, -1>> operator_series;

        bool wavepacket_finished;

        wavepacket(system_info<T, DIM>&, shell_input &);
        wavepacket();
        bool is_required();
        void printWavepacket();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();
};
