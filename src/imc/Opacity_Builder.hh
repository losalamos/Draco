//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.hh
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Opacity_Builder_hh__
#define __imctest_Opacity_Builder_hh__

//===========================================================================//
// class Opacity_Builder - 
//
// Purpose : build the Opacity and Mat_State objects, templated on MT
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Opacity.hh"
#include "imctest/Mat_State.hh"
#include "ds++/SP.hh"
#include <vector>

IMCSPACE

template<class MT>
class Opacity_Builder
{
private:
  // density, temperature, and opacity cell-centered fields
    typename MT::CCSF<double> rho;
    typename MT::CCSF<double> temp;
    typename MT::CCSF<double> opacity;

  // data received from XX_Interface
    vector<int> zone;
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> temperature;

  // Begin_Doc opacity_builder-int.tex
  // Begin_Verbatim 

public:
  // templated explicit constructor depends on interface type (IT)
    template<class IT>
    explicit Opacity_Builder(SP<IT>, SP<MT>);

  // build state member functions

  // build Mat_State helper functions
    SP< Mat_State<MT> > build_Mat();
    
  // build Opacity helper functions
    SP< Opacity<MT> > build_Opacity();	

  // End_Verbatim 
  // End_Doc 
};
    
CSPACE

#endif                          // __imctest_Opacity_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Opacity_Builder.hh
//---------------------------------------------------------------------------//
