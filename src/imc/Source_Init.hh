//----------------------------------*-C++-*----------------------------------//
// Source_Init.hh
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Source_Init_hh__
#define __imctest_Source_Init_hh__

//===========================================================================//
// class Source_Init - 
//
// Purpose : do source initialization for the source builder and the parallel 
//           plan
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Opacity.hh"
#include "imctest/Mat_State.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <string>

IMCSPACE

// draco components
using RNG::Rnd_Control;

// STL components
using std::string;

template<class MT>
class Source_Init
{
private:
  // data received from MT_Interface
    vector<double> evol_ext;
    vector<string> ss_pos;
    vector<double> ss_temp;
    vector<double> rad_temp;
    double delta_t;
    int npmax;
    double dnpdt;
    
  // source initialization data

  // number of particles for this cycle
    int npwant;

  // volume source variables
    typename MT::CCSF_double evol;
    double evoltot;

  // surface source variables
    typename MT::CCSF_double ess;
    typename MT::CCSF_int fss;
    double esstot;

  // radiation energy per cell, total
    typename MT::CCSF_double erad;
    double eradtot;

  // number of census particles per cell
    typename MT::CCSF_int ncen;
    int ncentot;

  // number of surface source and volume source particles
    typename MT::CCSF_int nvol;
    typename MT::CCSF_int nss;
    int nvoltot;
    int nsstot;

  // energy loss due to inadequate sampling of evol, ss, and initial census
    double eloss_vol;
    double eloss_ss;
    double eloss_cen;

  // energy weights for ss and vol emission source particles
    typename MT::CCSF_double ew_vol;
    typename MT::CCSF_double ew_ss;

  // maximum number of cells capable of fitting on a processor
    int capacity;

  // private member functions used to calc initial source information

  // number of source particles, census, source energies, number of volume
  // and surface sources
    void calc_initial_census(const MT &,const Opacity<MT> &, 
			     const Mat_State<MT> &, Rnd_Control &);
    void calc_source_energies(const Opacity<MT> &, const Mat_State<MT> &);
    void calc_source_numbers();

  // initial census service functions
    void calc_evol(const Opacity<MT> &, const Mat_State<MT> &);
    void calc_ess();
    void calc_erad();
    void calc_ncen_init();
    void write_initial_census(const MT &, Rnd_Control &);

public:
  // constructor
    template<class IT> Source_Init(SP<IT>, SP<MT>);

  // source initialyzer function
    void initialize(SP<MT>, SP<Opacity<MT> >, SP<Mat_State<MT> >, 
		    SP<Rnd_Control>, int);

  // accessor functions for Parallel_Builder
    int get_capacity() const { return capacity; }
    int get_ncentot() const { return ncentot; }
    int get_nsstot() const { return nsstot; }
    int get_nvoltot() const { return nvoltot; }
    int get_ncen(int cell) const { return ncen(cell); }
    int get_nvol(int cell) const { return nvol(cell); }
    int get_nss(int cell) const { return nss(cell); }
};

CSPACE

#endif                          // __imctest_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Source_Init.hh
//---------------------------------------------------------------------------//
