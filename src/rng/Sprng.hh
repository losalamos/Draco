//----------------------------------*-C++-*----------------------------------//
// Sprng.hh
// Thomas M. Evans
// Wed Apr 29 13:57:25 1998
//---------------------------------------------------------------------------//
// @> SPRNG random number class header file
//---------------------------------------------------------------------------//

#ifndef __rng_Sprng_hh__
#define __rng_Sprng_hh__

//===========================================================================//
// class Sprng - 
//
// Purpose : A random number class based on the SPRNG series of Generators
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "rng/Names.hh"
#include "ds++/Assert.hh"
#include "sprng.h"

RNGSPACE

class Sprng 
{
private:
    int *streamid;
    int streamnum;

public:
  // constructor
    Sprng(int *id, int snum) : streamid(id), streamnum(snum) {}
  // fake constructor for STL containers
    Sprng() { Check (0); } 

  // destructor, reclaim memory from SPRNG library
    ~Sprng() { free_sprng(streamid); }
   
  // get Random number
    double ran() const { return sprng(streamid); }

  // return the ID and number
    int* get_id() const { return streamid; }
    int get_num() const { return streamnum; }

  // do a diagnostic
    void print() const { print_sprng(streamid); }
};

CSPACE

#endif                          // __rng_Sprng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
