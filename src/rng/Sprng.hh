//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Sprng.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 29 13:57:25 1998
 * \brief  Sprng random number class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __rng_Sprng_hh__
#define __rng_Sprng_hh__

#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"

// Header file that includes the SPRNG library
#include <rng/config.h>
#include "rng_sprng.h"

#include <vector>
#include <cstdlib>

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
//
// revision history:
// -----------------
//  0) original
//  1)  4-30-98 : added reference counting to take care of memory in case the 
//                object is copied and deleted
//  2)  6-26-98 : added average tester 
//  3)  1-SEP-99: added doxygen comments
//  4) 20-DEC-01: added packing functions
// 
//===========================================================================//

class Sprng 
{
  private:
    // Reference Counting Class For Sprng Random Numbers.
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

    // Pointer to memory in Sprng library.
    SprngValue *streamid;

    // Number of this particular stream.
    int streamnum;

  public:
    // Constructors
    inline Sprng(int *, int);
    inline Sprng(const Sprng &);
    inline Sprng(const std::vector<char> &);

    // Destructor, reclaim memory from SPRNG library.
    inline ~Sprng();

    // Assignment operator.
    inline Sprng& operator=(const Sprng &);

    // Pack Sprng state.
    inline std::vector<char> pack() const;
   
    // >>> Services provided by Sprng class.

    // Get Random number.
    double ran() const { return sprng(streamid->id); }

    // Return the ID and number.
    int* get_id() const { return streamid->id; }
    int get_num() const { return streamnum; }

    // Test diagnostics.
    bool avg_test(int, double = .001) const;

    // Do a diagnostic.
    void print() const { print_sprng(streamid->id); }
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
    ++streamid->refcount;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 */
Sprng::Sprng(const std::vector<char> &packed)
    : streamid(0), streamnum(0)
{
    Require (packed.size() >= 2 * sizeof(int));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set the buffer
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the stream num and size of state
    int rng_size = 0;
    u >> streamnum >> rng_size;
    Check (streamnum >= 0);
    Check (rng_size >= 0);

    // unpack the random stream state
    char *prng = new char[rng_size];
    for (int i = 0; i < rng_size; i++)
	u >> prng[i];

    // now rebuild the sprng object
    int *rnid = unpack_sprng(prng);

    // now make a new streamid
    streamid = new SprngValue(rnid);

    // reclaim memory
    delete [] prng;

    Ensure (u.get_ptr() == &packed[0] + packed.size());
    Ensure (streamid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
Sprng::~Sprng()
{
    if (--streamid->refcount == 0)
	delete streamid;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a Sprng object into a vector<char>.
 */
std::vector<char> Sprng::pack() const
{
    Require (streamid);
    
    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number object and determine the size
    char *prng   = 0;
    int rng_size = pack_sprng(streamid->id, &prng);
    int size     = rng_size + 2 * sizeof(int);
    Check (prng);

    // now set the buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the stream number and rng size
    p << streamnum << rng_size;

    // pack the stream state
    for (int i = 0; i < rng_size; i++)
	p << prng[i];

    // free the prng buffer
    std::free(prng);

    Ensure (p.get_ptr() == &packed[0] + size);

    return packed;
}

//---------------------------------------------------------------------------//
// assignment operator

Sprng& Sprng::operator=(const Sprng &rhs)
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

} // end namespace rtt_rng

#endif                          // __rng_Sprng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
