//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Release.hh
 * \author Thomas M. Evans
 * \date   Thu May 27 15:24:01 1999
 * \brief  Release function for rng library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __rng_Release_hh__
#define __rng_Release_hh__

#include <string>

//===========================================================================//
/*!
 * \mainpage Overview of the RNG package
 * \version 1_4_0
 * \author Tom Evans, Todd Urbatsch
 *
 * The rng package is an object oriented random number generation package
 * that uses the SPRNG random number library.  
 */
//===========================================================================//
/*!
 * \namespace rtt_rng
 *
 * \brief Namespace that contains the rng package classes and variables.
 *
 * The rtt_rng namespace contains all classes and variables necessary to use
 * the rng package.  In particular, the namespace contains some "free"
 * variables that can be used to manage random number stream states.
 *
 * To use the rng package, all one has to do is include the Random.hh file in
 * their header and implementation files.  This file loads all rtt_rng
 * variables, the Sprng class, and the Rnd_Control class.  It is not
 * necessary to load these class headers individually.  
 */
//===========================================================================//

namespace rtt_rng 
{

//! Query package for the release number.
const std::string release();

} // end of rtt_rng

#endif                          // __rng_Release_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Release.hh
//---------------------------------------------------------------------------//
