//----------------------------------*-C++-*----------------------------------//
// Source_Init.hh
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Source_Init_hh__
#define __imc_Source_Init_hh__

//===========================================================================//
// class Source_Init - 
//
// Purpose : do source initialization for the source builder and the parallel 
//           plan
//
// revision history:
// -----------------
//  0) original
//  1)  7-13-98 : Added reproducible combing of census particles.
//  2)   8-3-98 : Temporary kluge (grrr) initializer for host code
//                interfacing 
//  3)   9-2-98 : Removed kluge (yeah!), fixed old_comb_census and
//                comb_census, the error was the additional particles 
//                (SP<PT> another) were being assigned to the original
//                particle.  Thus, these SP's were pointing to the same 
//                particle; the correct code is:
//                     SP<PT> another = new PT(*particle)
//  4) 6-18-99  : Removed the post-comb census adjustment.  Verified on 
//                five test problems.  Discussed in memo "Eliminating the 
//                Post-Comb Census Adjustment," XTM:99-49(U), July 13, 1999.
//  5) 10-6-99  : Added user-/host-defined surface source cells to calc_ess.
//  6) 10-27-99 : Added calculation and accessor functions for the fleck and
//                cummings time-explicit portion of the material volume
//                source. It will be used in the mat temp update in 
//                Global_Tally.
// 
//===========================================================================//

#include "Opacity.hh"
#include "Mat_State.hh"
#include "Particle.hh"
#include "Particle_Buffer.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <string>
#include <iostream>
#include <vector>

namespace rtt_imc 
{

// draco components
using rtt_rng::Rnd_Control;
using dsxx::SP;

// STL components
using std::string;
using std::vector;
using std::ostream;

template<class MT, class PT = Particle<MT> >
class Source_Init
{
private:
    // data received from MT_Interface
    vector<double> evol_ext;
    vector<double> rad_source;
    double rad_s_tend;
    vector<string> ss_pos;
    vector<double> ss_temp;
    vector< vector<int> > defined_surcells;
    vector<double> rad_temp;
    double delta_t;
    int npmax;
    int npnom;
    double dnpdt;
    string ss_dist;
    
    // source initialization data

    // number of particles for this cycle
    int npwant;

    // radiation volume emission source variables
    typename MT::CCSF_double evol;
    typename MT::CCSF_double evol_net;
    double evoltot;

    // external material volume source variables
    typename MT::CCSF_double mat_vol_src;
    double mat_vol_srctot;

    // surface source variables
    typename MT::CCSF_double ess;
    typename MT::CCSF_int fss;
    double esstot;

    // radiation energy per cell, total for census energy
    typename MT::CCSF_double ecen;
    double ecentot;

    // number of census particles per cell
    typename MT::CCSF_int ncen;
    int ncentot;
    SP<typename Particle_Buffer<PT>::Census> census;

    // number of surface source and volume source particles
    typename MT::CCSF_int nvol;
    typename MT::CCSF_int nss;
    int nvoltot;
    int nsstot;

    // energy loss due to inadequate sampling of evol, ss, and initial census
    double eloss_vol;
    double eloss_ss;
    double eloss_cen;

    // energy weights for census, ss, and vol emission source particles
    typename MT::CCSF_double ew_vol;
    typename MT::CCSF_double ew_ss;
    typename MT::CCSF_double ew_cen;

    // maximum number of cells capable of fitting on a processor
    int capacity;

    // slope of T_electron^4 in a cell -- using neighboring values
    typename MT::CCVF_double t4_slope;

    // private member functions used to calc initial source information

    // number of source particles, census, source energies, number of volume
    // and surface sources
    void calc_initial_census(const MT &, const Opacity<MT> &, 
			     const Mat_State<MT> &, Rnd_Control &, 
			     const int);
    void calc_source_energies(const Opacity<MT> &, const Mat_State<MT> &,
			      const int);
    void calc_source_numbers(const Opacity<MT> &, const int);
    void old_comb_census(const MT &, Rnd_Control &);
    void comb_census(const MT &, Rnd_Control &);

    // initial census service functions
    void calc_evol(const Opacity<MT> &, const Mat_State<MT> &, const int);
    void calc_ess();
    void calc_ecen();
    void calc_ncen_init();
    void write_initial_census(const MT &, Rnd_Control &);

    // calculate slope of T_electron^4 for volume emission
    void calc_t4_slope(const MT &, const Mat_State<MT> &);

public:
    // constructor
    template<class IT> Source_Init(SP<IT>, SP<MT>);

    // source initialyzer function
    void initialize(SP<MT>, SP<Opacity<MT> >, SP<Mat_State<MT> >, 
		    SP<Rnd_Control>, int);

    // accessor functions for Parallel_Builder
    int get_capacity() const { return capacity; }
    int get_nsstot() const { return nsstot; }
    int get_nvoltot() const { return nvoltot; }
    int get_nvol(int cell) const { return nvol(cell); }
    int get_nss(int cell) const { return nss(cell); }
    int get_fss(int cell) const { return fss(cell); }
    inline double get_t4_slope(int, int) const;
    double get_ew_vol(int cell) const { return ew_vol(cell); }
    double get_ew_ss(int cell) const { return ew_ss(cell); }
    string get_ss_dist() const { return ss_dist; }
    double get_evol_net(int cell) const { return evol_net(cell); }
    int num_cells() const { return ncen.get_Mesh().num_cells(); }

    // accessor functions for Global_Tally
    double get_delta_t() const { return delta_t; }
    double get_volume(int cell) const { return ncen.get_Mesh().volume(cell); }
    double get_mat_vol_src(int cell) const { return mat_vol_src(cell); }
    double get_mat_vol_srctot() const { return mat_vol_srctot; }

    // set and get functions for census stuff
    void set_ncen(int cell, int num) { ncen(cell) = num; }
    void set_ecen(int cell, double cen) { ecen(cell) = cen; }
    int get_ncen(int cell) const { return ncen(cell); }
    double get_ecen(int cell) const { return ecen(cell); }
    void set_ecentot(double cent) { ecentot = cent; }
    void set_ncentot(int num) { ncentot = num; }
    int get_ncentot() const { return ncentot; }
    inline void set_census(SP<typename Particle_Buffer<PT>::Census>);
    inline SP<typename Particle_Buffer<PT>::Census> get_census() const;	

    // get functions for energy
    double get_evoltot() const { return evoltot; }
    double get_esstot() const { return esstot; } 
    double get_ecentot() const { return ecentot; }
    double get_eloss_ss() const { return eloss_ss; }
    double get_eloss_vol() const { return eloss_vol; }
    double get_eloss_cen() const { return eloss_cen; }

    // diagnostic functions
    void print(ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT, class PT>
inline ostream& operator<<(ostream &out, 
			   const Source_Init<MT,PT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// inline functions for source_init
//---------------------------------------------------------------------------//
// set the census for updates between cycles

template<class MT, class PT> inline 
void Source_Init<MT,PT>::set_census(SP<typename Particle_Buffer<PT>::Census> 
				    census_)
{
    // we must update this with a valid census
    Require (census_);
    census = census_;
}

//---------------------------------------------------------------------------//
// return the census

template<class MT, class PT> inline
SP<typename Particle_Buffer<PT>::Census> Source_Init<MT,PT>::get_census()
    const 
{
    return census;
}

//---------------------------------------------------------------------------//
// return the slope of the temperature to the 4th power in each cell

template<class MT, class PT>
inline double Source_Init<MT,PT>::get_t4_slope(int dim, int cell) const 
{ 
    return t4_slope(dim, cell); 
}

} // end namespace rtt_imc

#endif                          // __imc_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Init.hh
//---------------------------------------------------------------------------//
