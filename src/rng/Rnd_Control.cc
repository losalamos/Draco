//----------------------------------*-C++-*----------------------------------//
// Rnd_Control.cc
// Thomas M. Evans
// Wed Apr 29 16:08:59 1998
//---------------------------------------------------------------------------//
// @> Rnd_Control class implementation file
//---------------------------------------------------------------------------//

#include "rng/Rnd_Control.hh"
#include "ds++/Assert.hh"

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
// spawn a new random number stream

SP<Sprng> Rnd_Control::spawn(Sprng &random)
{
  // declare variables necessary to spawn a stream
    int **newstream;
    int numspawn;
    
  // spawn a new stream
    numspawn = spawn_sprng(random.get_id(), 1, &newstream);
    Check (numspawn == 1);

  // create a new SPRNG random number object with new stream
    SP<Sprng> ran = new Sprng(newstream[0], random.get_num());

  // retrun the new random object
    return ran;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Rnd_Control.cc
//---------------------------------------------------------------------------//
