//----------------------------------*-C++-*----------------------------------//
// Constants.hh
// Thomas M. Evans
// Fri Jan 30 16:23:53 1998
//---------------------------------------------------------------------------//
// @> IMC::Global namespace constants
//---------------------------------------------------------------------------//

#ifndef __imc_Constants_hh__
#define __imc_Constants_hh__

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
// 
//===========================================================================//

#include "imc/Names.hh"
#include <cfloat>

IMCSPACE
GLOBALSPACE

//---------------------------------------------------------------------------//
// FUNDAMENTAL CONSTANTS
//---------------------------------------------------------------------------//

const double pi = 3.14159265359;
const double huge = DBL_MAX;
const double epsilon = DBL_EPSILON;

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

CSPACE
CSPACE

#endif                          // __imc_Constants_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Constants.hh
//---------------------------------------------------------------------------//
