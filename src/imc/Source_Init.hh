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
#include "ds++/SP.hh"
#include <string>

IMCSPACE

using std::string;

template<class MT>
class Source_Init
{
private:
  // data received from MT_Interface
    vector<double> evol_ext;
    vector<string> ss_pos;
    vector<double> ss_temp;
    double delta_t
    
  // source initialization data

  // volume source variables
    typename MT::CCSF<double> evol;
    double evoltot;

  // surface source variables
    typename MT::CCVF<double> ess;
    double esstot;


  // private member functions used to calc initial source information

  // number of source particles, census, source energies, number of volume
  // and surface sources
    void Calc_num_part();
    void Calc_initial_census(const Opacity<MT> &, const Mat_State<MT> &);
    void Calc_source_energies();
    void Calc_source_numbers();

  // initial census service functions
    void Calc_evol(const Opacity<MT> &, const Mat_State<MT> &);
    void Calc_ess();
    void Calc_erad();
    void Calc_ncen_init();
    void Write_initial_census();

public:
    template<class IT>
    explicit Source_Init(SP<IT>, SP<MT>);

  // source initialyzer function
    void Initialize(const Opacity<MT> &, const Mat_State<MT> &);

CSPACE

#endif                          // __imctest_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Source_Init.hh
//---------------------------------------------------------------------------//
