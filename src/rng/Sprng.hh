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
// Purpose : A random number class based on the SPRNG series of Generators.
//
// revision history:
// -----------------
//  0) original
//  1)  4-30-98 : added reference counting to take care of memory in case the 
//                object is copied and deleted
//  2)  6-26-98 : added average tester 
// 
//===========================================================================//

#include "Names.hh"
#include "ds++/Assert.hh"

// Header file that includes the SPRNG library
#include <rng/config.h>
#include "rng_sprng.h"

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
    inline Sprng();

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

  // test diagnostics
    bool avg_test(int, double = .001) const;

  // do a diagnostic
    void print() const { print_sprng(streamid->id); }
};

//---------------------------------------------------------------------------//
// inline members for Sprng
//---------------------------------------------------------------------------//
// constructor

inline Sprng::Sprng(int *idval, int number)
    : streamid(new SprngValue(idval)), streamnum(number) {}

//---------------------------------------------------------------------------//
// copy constructor

inline Sprng::Sprng(const Sprng &rhs)
    : streamid(rhs.streamid), streamnum(rhs.streamnum)
{
    ++streamid->refcount;
}

//---------------------------------------------------------------------------//
// destructor

inline Sprng::~Sprng()
{
    if (--streamid->refcount == 0)
	delete streamid;
}

//---------------------------------------------------------------------------//
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
    streamid  = rhs.streamid;
    streamnum = rhs.streamnum; 
    ++streamid->refcount;
    return *this;
}

//---------------------------------------------------------------------------//
// default constructor for STL classes

inline Sprng::Sprng() 
    : streamid(0) 
{ 
  // constructor for use with STL containers, this cannot be used
    Insist (0, "You tried to default construct a Sprng!"); 
} 

CSPACE

#endif                          // __rng_Sprng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
