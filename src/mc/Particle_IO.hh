//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_IO.hh
 * \author Thomas M. Evans
 * \date   Fri Dec 21 10:10:55 2001
 * \brief  Particle_IO class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Particle_IO_hh__
#define __mc_Particle_IO_hh__

#include "Particle_Stack.hh"
#include <iostream>

namespace rtt_mc
{

// Forward declarations of Particle_Buffer.
template<class PT> class Particle_Buffer;
 
//===========================================================================//
/*!
 * \class Particle_IO
 *
 * Do Particle_IO for file manipulation.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class PT>
class Particle_IO 
{
  public:
    // Useful typedefs.
    typedef std::ostream                           std_ostream;
    typedef std::istream                           std_istream;
    typedef typename Particle_Containers<PT>::Bank Bank;

  public:

    // >>> I/O FUNCTIONS.

    // Write a particle to disk.
    static void write_particle(std_ostream &, const PT &);

    // Write a Particle_Buffer to disk.
    static void write_Particle_Buffer(std_ostream &, 
				      const Particle_Buffer<PT> &);

    // Read all particles on the file
    static void read_particles(std_istream &, Bank &);
};

} // end namespace rtt_mc

#endif                          // __mc_Particle_IO_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Particle_IO.hh
//---------------------------------------------------------------------------//
