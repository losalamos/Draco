//----------------------------------*-C++-*----------------------------------//
// Rnd_Control.cc
// Thomas M. Evans
// Wed Apr 29 16:08:59 1998
//---------------------------------------------------------------------------//
// @> Rnd_Control class implementation file
//---------------------------------------------------------------------------//

#include "rng/Rnd_Control.hh"

RNGSPACE

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//
// get a new Random number object with its own stream

SP<Sprng> Rnd_Control::get_rn()
{
  // declare a stream
    int *id = init_sprng(streamnum, number, seed, parameter);

  // create a new Rnd object
    SP<Sprng> random = new Sprng(id, number);

  // advance the counter
    streamnum++;

  // return the object
    return random;
}

//---------------------------------------------------------------------------//
// spawn a new random

//---------------------------------------------------------------------------//
//                              end of Rnd_Control.cc
//---------------------------------------------------------------------------//
