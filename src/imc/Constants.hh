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
//===========================================================================//

#include "imc/Names.hh"
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

//---------------------------------------------------------------------------//
// MAXIMUM BUFFER SIZES
//---------------------------------------------------------------------------//
// maximum number of Particles in a stack
const int buffer_s = 1000;

// maximum size of Particle double stack
const int buffer_d = buffer_s * (3 + 3 + 1 + 1 + 1);

// maximum size of Particle int stack
const int buffer_i = buffer_s * (1 + 1);

// maximum size of Particle char stack
const int buffer_c = buffer_s * (500);

//---------------------------------------------------------------------------//
// Physical constants
//---------------------------------------------------------------------------//

// radiation constant
const double a = 0.01372;

// speed of light (m/s)
const double c = 2.99792458e8;


CSPACE
CSPACE

#endif                          // __imc_Constants_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Constants.hh
//---------------------------------------------------------------------------//
