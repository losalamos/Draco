//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.hh
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Opacity_Builder_hh__
#define __imc_Opacity_Builder_hh__

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

#include "imc/Names.hh"
#include "imc/Opacity.hh"
#include "imc/Mat_State.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

IMCSPACE

using std::string;

template<class MT>
class Opacity_Builder
{
private:

  // data received from XX_Interface
    vector<int> zone;
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> temperature;
    vector<double> specific_heat;
    double implicitness;
    double delta_t;
    string analytic_opacity;

  // Begin_Doc opacity_builder-int.tex
  // Begin_Verbatim 

public:
  // templated explicit constructor depends on interface type (IT)
    template<class IT>
    explicit Opacity_Builder(SP<IT>);

  // build state member functions

  // build Mat_State helper functions
    SP< Mat_State<MT> > build_Mat(SP<MT>);
    
  // build Opacity helper functions
    SP< Opacity<MT> > build_Opacity(SP<MT>, SP<Mat_State<MT> >);	

  // End_Verbatim 
  // End_Doc 
};
    
CSPACE

#endif                          // __imc_Opacity_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity_Builder.hh
//---------------------------------------------------------------------------//
