//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Sprng.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 29 13:57:25 1998
 * \brief  Sprng random number class header file.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_rng_Sprng_hh
#define rtt_rng_Sprng_hh

#include <vector>
#include "ds++/Assert.hh"

// Header file that includes the SPRNG library
#include "rng_sprng.h"

namespace rtt_rng 
{

//===========================================================================//
/*!
 * \class Sprng
 *
 * \brief A random object class that utilizes the SPRNG random number
 *        library.
 * 
 * The Sprng random number class is a wrapper class for the <a
 * href="http://www.ncsa.uiuc.edu/Apps/SPRNG/">SPRNG (Scalable Parallel
 * Random Number Generator)</a>) libraries.  The class is nominally
 * controlled by the Rnd_Control class; however, Sprng random number objects
 * can be created directly from native SPRNG functions.  See the \link
 * tstSprng.cc Sprng examples \endlink for more detail.
 *
 * The primary purpose of the Sprng class is to manage the memory that is
 * controlled by the SPRNG library. Sprng accomplishes this through reference
 * counting.  Sprng objects can be copied an assigned freely without worrying
 * about memory.  When the last Sprng object that accesses a particular
 * random number state is destroyed, the memory controlling that state is
 * released.  Note that when a Sprng is copied, both instantiations of the
 * Sprng object access the same random number state.  See the \link
 * tstSprng.cc examples \endlink for more info.
 *
 * \sa <a href="http://www.ncsa.uiuc.edu/Apps/SPRNG/">SPRNG</a>,
 * rtt_rng, Rnd_Control, Random.hh 
 */
/*!
 * \example rng/test/tstSprng.cc
 *
 * Example of Sprng class usage.  The Sprng random numbers are created by
 * direct usage of the SPRNG random number library functions.  A more
 * preferred way of creating Sprng objects is to use the Rnd_Control class
 * that "hides" SPRNG library details from the user.  
 */
//===========================================================================//

class Sprng 
{
  private:
    // Reference Counting Class For Sprng Random Numbers.
    struct SprngValue
    {
	// counter and id to Sprng library
	int *id;
	int refcount;

	// constructor
	SprngValue(int *idnum) :  id(idnum), refcount(1) {}

	// destructor
	~SprngValue() { free_sprng(id); }
    };	

    // Pointer to memory in Sprng library.
    SprngValue *streamid;

    // Number of this particular stream.
    int streamnum;

    // Size of the packed state
    static int packed_size;

  public:
    // Constructors
    inline Sprng() : streamid(0) {}
    inline Sprng(int *, int);
    inline Sprng(const Sprng &);
    Sprng(const std::vector<char> &);

    // Destructor, reclaim memory from SPRNG library.
    inline ~Sprng();

    // Assignment operator.
    Sprng& operator=(const Sprng &);

    // Pack Sprng state.
    std::vector<char> pack() const;
    void pack(std::vector<char>& packed, std::vector<int>& tmpbuf) const;
   

    // >>> Services provided by Sprng class.

    // Get Random number.
    double ran() const { Require(streamid); return sprng(streamid->id); }

    // Return the ID and number.
    int* get_id() const { Require(streamid); return streamid->id; }
    int get_num() const { Require(streamid); return streamnum; }

    // Return the packed size
    int get_size() const;

    // Test diagnostics.
    bool avg_test(int, double = .001) const;

    // Do a diagnostic.
    void print() const { Require(streamid); print_sprng(streamid->id); }
};

//---------------------------------------------------------------------------//
// INLINE SPRNG MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Sprng::Sprng(int *idval, int number)
    : streamid(new SprngValue(idval)), streamnum(number) {}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
Sprng::Sprng(const Sprng &rhs)
    : streamid(rhs.streamid), streamnum(rhs.streamnum)
{
    if(streamid)
        ++streamid->refcount;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
Sprng::~Sprng()
{
    if (streamid && --streamid->refcount == 0)
	delete streamid;
}

} // end namespace rtt_rng

#endif                          // rtt_rng_Sprng_hh

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
