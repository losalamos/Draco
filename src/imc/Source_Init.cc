//----------------------------------*-C++-*----------------------------------//
// Source_Init.cc
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Source_Init.hh"

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// source initialyzer
template<class MT>
void Source_Init<MT>::Initialize();
{
  // calculate number of particles
    Calc_np();

  // on first pass do initial source census
    if (cycle == 1)
	Calc_initial_census();

  // calculate source energies
    Calc_energies();

  // calculate source numbers
    Calc_numbers();
}

//---------------------------------------------------------------------------//
//                              end of Source_Init.cc
//---------------------------------------------------------------------------//
