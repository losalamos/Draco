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

namespace rtt_rng {
class Rnd_Control;
}

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
    int *sid;

    // Size of the packed state
    static int packed_size;

  public:
    // Constructors
    inline Sprng(int *, int);
    inline Sprng(const Sprng &);
    Sprng(const std::vector<char> &);

    // Destructor, reclaim memory from SPRNG library.
    inline ~Sprng();

    // Assignment operator.
    Sprng& operator=(const Sprng &);

    // Pack Sprng state.
    std::vector<char> pack() const;
   
    // >>> Services provided by Sprng class.

    // Get Random number.
    double ran() const { return sprng(sid); }

    // Return the ID and number.
    int* get_id() const { return sid; }
    int get_num() const { return streamnum; }

    // Return the packed size
    int get_size() const;

    // Test diagnostics.
    bool avg_test(int, double = .001) const;

    // Do a diagnostic.
    void print() const { print_sprng(sid); }
};

//---------------------------------------------------------------------------//
// INLINE SPRNG MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Sprng::Sprng(int *idval, int number)
    : streamid(new SprngValue(idval)), streamnum(number), sid(streamid->id) {}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
Sprng::Sprng(const Sprng &rhs)
    : streamid(rhs.streamid), streamnum(rhs.streamnum), sid(streamid->id) 
{
    ++streamid->refcount;
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



} // end namespace rtt_rng

#endif                          // __rng_Sprng_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Sprng.hh
//---------------------------------------------------------------------------//
