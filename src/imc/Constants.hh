//----------------------------------*-C++-*----------------------------------//
// Constants.hh
// Thomas M. Evans
// Fri Jan 30 16:23:53 1998
//---------------------------------------------------------------------------//
// @> IMC::Global namespace constants
//---------------------------------------------------------------------------//

#ifndef __imctest_Constants_hh__
#define __imctest_Constants_hh__

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
//===========================================================================//

#include "imctest/Names.hh"
#include <cmath>
#include <cfloat>

IMCSPACE
GLOBALSPACE

//---------------------------------------------------------------------------//
// FUNDAMENTAL CONSTANTS
//---------------------------------------------------------------------------//
const double pi = 2.0 * std::asin(1.0);
const double huge = DBL_MAX;
const double epsilon = DBL_EPSILON;

CSPACE
CSPACE

#endif                          // __imctest_Constants_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Constants.hh
//---------------------------------------------------------------------------//
