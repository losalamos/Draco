//----------------------------------*-C++-*----------------------------------//
// Rnd_Control.cc
// Thomas M. Evans
// Wed Apr 29 16:08:59 1998
//---------------------------------------------------------------------------//
// @> Rnd_Control class implementation file
//---------------------------------------------------------------------------//

#include "Rnd_Control.hh"
#include "ds++/Assert.hh"
#include <cstdlib>

// header file for SPRNG package
#include <rng/config.h>
#include "rng_sprng.h"

RNGSPACE

using std::free;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
// calculates the size of a stored Random Number state

Rnd_Control::Rnd_Control(int s, int n, int sn, int p)
    : seed(s), number(n), streamnum(sn), parameter(p)
{
  // calculate the size of a stored random number state
    int  *id;
    char *buffer;

  // operations
    id   = init_sprng(0, number, seed, parameter);
    size = pack_sprng(id, &buffer);
    
  // deallocate memory
    free(buffer);
    free_sprng(id);
}

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//
// get a new Random number object with its own stream

Sprng Rnd_Control::get_rn()
{
  // declare a stream
    int *id = init_sprng(streamnum, number, seed, parameter);

  // create a new Rnd object
    Sprng random(id, streamnum);

  // advance the counter
    streamnum++;

  // return the object
    return random;
}

//---------------------------------------------------------------------------//
// get a new Random number object and reset the stream

Sprng Rnd_Control::get_rn(int snum)
{
  // reset streamnum
    streamnum = snum;

  // declare a stream
    int *id = init_sprng(streamnum, number, seed, parameter);

  // create a new Rnd object
    Sprng random(id, streamnum);

  // advance the counter
    streamnum++;

  // return the object
    return random;
}

//---------------------------------------------------------------------------//
// spawn a new random number stream

Sprng Rnd_Control::spawn(const Sprng &random) const
{
  // declare variables necessary to spawn a stream
    int **newstream;
    int numspawn;
    
  // spawn a new stream
    numspawn = spawn_sprng(random.get_id(), 1, &newstream);
    Check (numspawn == 1);

  // create a new SPRNG random number object with new stream
    Sprng ran(newstream[0], random.get_num());

  // free some memory
    free(newstream);

  // return the new random object
    return ran;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Rnd_Control.cc
//---------------------------------------------------------------------------//
