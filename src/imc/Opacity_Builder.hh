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

#include "Names.hh"
#include "OS_Parser.hh"
#include "SP.hh"
#include <vector>

IMCSPACE

template<class MT>
class Opacity_Builder
{
private:
  // data received from XX_Parser
    vector<int> zone;
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> temperature;

public:
  // templated explicit constructor depends on parser type (PT)
    template<class PT>
    explicit Opacity_Builder(SP<PT> parser)
    {
      // assign data members from the parser
	zone        = parser->Zone();
	mat_zone    = parser->Mat_zone();
	density     = parser->Density();
	kappa       = parser->Kappa();
	temperature = parser->Temperature();
    }

  // build state member functions

  // build Mat_State helper functions
    void Build_Mat(SP<MT>);
    
  // build Opacity helper functions
    void Build_Opacity(SP<MT>);	
};
    
CSPACE

#endif                          // __imctest_Opacity_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Opacity_Builder.hh
//---------------------------------------------------------------------------//
