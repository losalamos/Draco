//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Rnd_Control.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 29 16:08:59 1998
 * \brief  Rnd_Control class implementation file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <cstdlib>
#include <vector>
#include "Rnd_Control.hh"

namespace rtt_rng 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!  
 * \brief Rnd_Control class constructor.
 *
 * Constructor for the Rnd_Control class only requires one argument with 3
 * additional options.
 *
 * \param s seed value 
 * \param n total number of independent streams; default = 1.0e9 
 * \param sn stream index; default = 0 
 * \param p Sprng parameter that determines the size of the Sprng state, see
 *        the Sprng library documentation for details; default = 1 
 */
Rnd_Control::Rnd_Control(int s, int n, int sn, int p)
    : seed(s), number(n), streamnum(sn), parameter(p)
{
    Require (streamnum <= number);

    // make a spring object and pack it to determine the size
    Sprng temp = make_random_number_generator();

    // pack it and set size
    std::vector<char> pack = temp.pack();
    size                   = pack.size();
    Check (size >= 0);
    Check (streamnum == sn);
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Create a Sprng random number object at the current stream index.
 * 
 * \return  Sprng random number object
 *
 * This function creates a Sprng random number state at the current stream
 * index.  After creation, the stream index is automatically advanced by one.
 * Additionally, the created Sprng object takes ownership of the random
 * number state from the Sprng library.  No further work on the part of
 * Rnd_Control is required.  
 */
Sprng Rnd_Control::get_rn()
{
    Require (streamnum <= number);

    // create a new Rnd object
    Sprng random = make_random_number_generator();

    // advance the counter
    streamnum++;

    // return the object
    return random;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a Sprng random number object at a given stream index.
 *
 * \param snum a user-requested random number index
 * \return Sprng random number object
 *
 * This function is identical to rtt_rng::Rnd_Control::get_rn() with the
 * exception: the user enters a stream index for the random number state.
 * Note that the entered stream index becomes the new Rnd_Control stream
 * index from which future random numbers are created.  The entered stream
 * index is incremented by one after creating the random number. 
 */
Sprng Rnd_Control::get_rn(int snum)
{
    Require (snum <= number);

    // reset streamnum
    streamnum = snum;

    // create a new Rnd object
    Sprng random = make_random_number_generator();

    // advance the counter
    streamnum++;

    // return the object
    return random;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Spawn a new random number object.
 *
 * \param random valid Sprng random number object
 * \return Sprng random number object
 *
 * This function spawns a new, independent random number object from an
 * existing Sprng random number object.  The new Sprng random number state
 * has as its index the index of the parent Sprng object from which it has
 * been spawned.  
 * 
 * Memory control of the spawned stream is given to the returned Sprng
 * object.  The spawned object is--from this point on--independent of the
 * stream from which it spawned.  This class is used to guarantee
 * reproducibility in random number codes.  By spawning only from existing
 * random number states, reproducibility in numeric operations can be
 * achieved.  
*/
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
    std::free(newstream);

    // return the new random object
    return ran;
}

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                              end of Rnd_Control.cc
//---------------------------------------------------------------------------//
