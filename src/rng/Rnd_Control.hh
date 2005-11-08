//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Rnd_Control.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 29 16:08:58 1998
 * \brief  Rnd_Control header file.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_rng_Rnd_Control_hh
#define rtt_rng_Rnd_Control_hh

// header file for SPRNG package
#include "rng_sprng.h"

#include "ds++/Assert.hh"
#include "Sprng.hh"

namespace rtt_rng 
{

//===========================================================================//
/*!
 * \class Rnd_Control 
 *
 * \brief Initialize and control the creation of Sprng random rumber objects.
 * 
 * The Rnd_Control class is used to manage the Sprng random number classes.
 * The class' main function is to create a new Sprng random number object.
 * By using the Rnd_Control class, users are not required to "ever" access
 * the SPRNG random number library explicitly. Rnd_Control calls all the
 * appropriate SPRNG library functions to set up a random number state.
 * Ownership of the random number state is maintained by the created Sprng
 * class object. By using the controller, all memory management \i vis-a-vis
 * the SPRNG library is taken care of automatically.
 *
 * The usage is simple, make a Rnd_Control object and then query it for Sprng
 * random number objects.  Once a Sprng object is built, memory management of
 * the Sprng library associated with that state is taken care of by the Sprng
 * object.  The Rnd_Control object keeps track of the number of independent *
 * streams created through the private member Rnd_Control::streamnum (stream
 * * index).  This index is automatically incremented by one everytime a
 * random number is created. The linear progression of random number states
 * can be interupted by giving an optional stream index (integer) to the
 * get_rn(int) function or by reseting the random number stream index through
 * the set_num(int) function.
 *
 * \sa <a href="http://www.ncsa.uiuc.edu/Apps/SPRNG/">SPRNG (Scalable
 * Parallel Random Number Generator Library)</a>, Sprng, Random.hh, rtt_rng 
 */
//===========================================================================//

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

  private:
    // >>> Implementation.
    inline Sprng make_random_number_generator() const;

  public:
    // Constructor.
    Rnd_Control(int s, int n = 1000000000, int sn = 0, int p = 1);

    // Create Sprng objects.
    Sprng get_rn();
    Sprng get_rn(int);

    // Spawn a new random number object.
    Sprng spawn(const Sprng &) const;

    //! Query for the current random number stream index.
    int get_num() const { return streamnum; }

    //! Set (reset) the random number stream index.
    void set_num(int num) { streamnum = num; }

    //! Query size of a packed random number state.
    int get_size() const { return size; }

    //! Get the seed value used to initialize the SPRNG library.
    int get_seed() const { return seed; }

    //! Return the total number of current streams set.
    int get_number() const { return number; }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Make a Sprng random number object.
 *
 * This function is used by the public interface to get a Sprng random number
 * object.
 */
Sprng Rnd_Control::make_random_number_generator() const 
{
    Require (streamnum <= number);

    // declare a stream
    int *id = init_sprng(streamnum, number, seed, parameter);

    // create a new Rnd object
    Sprng random(id, streamnum);

    // return the object
    return random; 
}

} // end namespace rtt_rng

#endif                          // rtt_rng_Rnd_Control_hh

//---------------------------------------------------------------------------//
//                              end of rng/Rnd_Control.hh
//---------------------------------------------------------------------------//
