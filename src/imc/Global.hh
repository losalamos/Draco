//----------------------------------*-C++-*----------------------------------//
// Global.hh
// Thomas M. Evans
// Mon Jun  8 11:29:48 1998
//---------------------------------------------------------------------------//
// @> Global namespace declarations
//---------------------------------------------------------------------------//

#ifndef __imc_Global_hh__
#define __imc_Global_hh__

//===========================================================================//
// namespace Global - 
//
// Purpose : groups all functions, variables, and constants that are defined 
//           in the IMC::Global namespace into one header.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Math.hh"
#include "imc/Constants.hh"

IMCSPACE
GLOBALSPACE

//---------------------------------------------------------------------------//
// PROBLEM VARIABLES
//---------------------------------------------------------------------------//

extern int rn_stream;

//---------------------------------------------------------------------------//
// MAXIMUM BUFFER SIZES
//---------------------------------------------------------------------------//

// max number of double data
const int data_d = 3 + 3 + 1 + 1 + 1;

// max number of int data
const int data_i = 1 + 1;

// maximum number of Particles in a stack
extern int buffer_s;
const int default_s = 1000;

// maximum size of Particle double stack
extern int buffer_d;

// maximum size of Particle int stack
extern int buffer_i;

// maximum size of Particle char stack
extern int buffer_c;

CSPACE
CSPACE

#endif                          // __imc_Global_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Global.hh
//---------------------------------------------------------------------------//
