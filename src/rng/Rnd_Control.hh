//----------------------------------*-C++-*----------------------------------//
// Rnd_Control.hh
// Thomas M. Evans
// Wed Apr 29 16:08:58 1998
//---------------------------------------------------------------------------//
// @> Rnd_Control header file
//---------------------------------------------------------------------------//

#ifndef __rng_Rnd_Control_hh__
#define __rng_Rnd_Control_hh__

//===========================================================================//
// class Rnd_Control - 
//
// Purpose : initialize and control the creation of SPRNG Random Number
//           objects
//
// revision history:
// -----------------
//  0) original
//  1)  5-21-98 : added get_size() function that determines the size of the
//                stored random number state
// 
//===========================================================================//

#include "Names.hh"
#include "Sprng.hh"

RNGSPACE

class Rnd_Control 
{
private:
  // seed for initialization of random number streams
    int seed;
  // total number of streams
    int number;
  // number of current stream
    int streamnum;
  // control parameter for stream inits
    int parameter;
  // size of packed stream state
    int size;

public:
  // constructor
    Rnd_Control(int, int = 1000000000, int = 0, int = 1);

  // start a new random number stream
    Sprng get_rn();
    Sprng get_rn(int);

  // spawn a new random number stream
    Sprng spawn(const Sprng &) const;

  // query for the number of random streams
    int get_num() const { return streamnum; }
    void set_num(int num) { streamnum = num; }

  // query size of packed random number state
    int get_size() const { return size; }

  // get the seed value and total number of current streams set here
    int get_seed() const { return seed; }
    int get_number() const { return number; }
};

CSPACE

#endif                          // __rng_Rnd_Control_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Rnd_Control.hh
//---------------------------------------------------------------------------//
