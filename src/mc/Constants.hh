//----------------------------------*-C++-*----------------------------------//
// Constants.hh
// Thomas M. Evans
// Fri Jan 30 16:23:53 1998
//---------------------------------------------------------------------------//
// @> mc::global namespace constants
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Constants_hh
#define rtt_mc_Constants_hh

//===========================================================================//
// Namespace Constants - 
//
// Purpose : provide constants which are used globally
//
// revision history:
// -----------------
//  0) original
//  1)   2-4-98 : added huge and epsilon constants
//  2)  3-20-98 : split former Global class into Constants and Math
//                namespaces, moved Global namespace inside of IMC namespace
//
//  3)   6-5-98 : added variables that are needed Globally
//  4)   6-8-98 : reorganized so that all Global stuff is called from
//                Global.hh and the variables are placed in Global.cc
//  5)  4-13-99 : moved to mc package
// 
//===========================================================================//

#include <cmath>
#include <cfloat>

namespace rtt_mc 
{
namespace global
{

//---------------------------------------------------------------------------//
// FUNDAMENTAL CONSTANTS
//---------------------------------------------------------------------------//

const double pi       = 2.0 * std::asin(1.0);
const double huge     = 1.0e22;
const int    huge_int = 2000000000;
const double epsilon  = DBL_EPSILON;

//---------------------------------------------------------------------------//
// Physical constants
//---------------------------------------------------------------------------//

// radiation constant, jerks/cm^3/keV^4
const double a = 0.01372;

// speed of light, cm/shake
const double c = 2.99792458e2;

// boltzman constant, keV/K
const double k = 8.617308e-8;

// Planck constant, keV-s
const double h = 4.135706e-18;

} // end namespace global
} // end namespace rtt_mc

#endif                          // rtt_mc_Constants_hh

//---------------------------------------------------------------------------//
//                              end of mc/Constants.hh
//---------------------------------------------------------------------------//
