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
// 0) original
// 
//===========================================================================//

#include "rng/Names.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "sprng.h"

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

public:
  // constructor
    Rnd_Control(int s, int n = 1000000, int sn = 0, int p = 1) 
	: seed(s), number(n), streamnum(sn), parameter(p) {}

  // start a new random number stream
    Sprng get_rn();
    Sprng get_rn(int);

  // spawn a new random number stream
    Sprng spawn(Sprng &);

  // query for the number of random streams
    int get_num() const { return streamnum; }
    void set_num(int num) { streamnum = num; }
};

CSPACE

#endif                          // __rng_Rnd_Control_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Rnd_Control.hh
//---------------------------------------------------------------------------//
