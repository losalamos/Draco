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
//  0)  original
//  1)  4-30-98: added reference counting to take care of memory in case the 
//               object is copied and deleted
// 
//===========================================================================//

#include "rng/Names.hh"
#include "ds++/Assert.hh"
#include "sprng.h"

RNGSPACE

class Sprng 
{
private:
  // reference counting class for SPRNG random numbers
    struct SprngValue
    {
      // counter and id to Sprng library
	int refcount;
	int *id;

      // constructor
	SprngValue(int *idnum) : refcount(1), id(idnum) {}

      // destructor
	~SprngValue() { free_sprng(id); }
    };	

  // pointer to memory in Sprng library
    SprngValue *streamid;
  // number of this particular stream
    int streamnum;

public:
  // creation, copy, and assignment

  // constructors
    inline Sprng(int *, int);
    inline Sprng(const Sprng &);

  // fake constructor for STL containers
    Sprng() : streamid(0) { Check (0); } 

  // destructor, reclaim memory from SPRNG library
    inline ~Sprng();

  // assignment operator
    inline Sprng& operator=(const Sprng &);
   
  // services provided by Sprng class

  // get Random number
    double ran() const { return sprng(streamid->id); }

  // return the ID and number
    int* get_id() const { return streamid->id; }
    int get_num() const { return streamnum; }

  // do a diagnostic
    void print() const { print_sprng(streamid->id); }
};

//---------------------------------------------------------------------------//
// inline members for Sprng
//---------------------------------------------------------------------------//

// constructors
inline Sprng::Sprng(int *idval, int number)
    : streamid(new SprngValue(idval)), streamnum(number) {}

inline Sprng::Sprng(const Sprng &rhs)
    : streamid(rhs.streamid), streamnum(rhs.streamnum)
{
    ++streamid->refcount;
}

// destructor
inline Sprng::~Sprng()
{
    if (--streamid->refcount == 0)
	delete streamid;
}

// assignment operator
inline Sprng& Sprng::operator=(const Sprng &rhs)
{
  // check to see if the values are the same
    if (streamid == rhs.streamid && streamnum == rhs.streamnum)
	return *this;

  // destroy this' value if it was the last random number in a particular
  // stream 
    if (--streamid->refcount == 0)
	delete streamid;

  // do assignment
    streamid = rhs.streamid;
    ++streamid->refcount;
    return *this;
}

CSPACE

#endif                          // __rng_Sprng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
