//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Release.hh
 * \author Thomas M. Evans
 * \date   Tue Apr 13 17:58:54 1999
 * \brief  Release function for mc library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Release_hh__
#define __mc_Release_hh__

//===========================================================================//
/*!
 * \page mc_overview Overview of the MC package
 * \version 1_2_0
 * \author Tom Evans, Todd Urbatsch
 *
 * The mc package is a set of general Monte Carlo (mc) components that can be
 * used in all X-6 Monte Carlo codes.  In particular it includes the
 * following types of classes:
 *
 * \arg \b Meshes: mesh types (MT) that are designed for Monte Carlo use \arg
 * \b Topology: Topology class that is used for spatial mesh decompositions
 * \arg \b Particle Base: base particle classes for Monte Carlo transport
 * \arg \b Particle_Buffer: particle containers for storage and communication
 * across processors.  
 */
//===========================================================================//
/*!
 * \namespace rtt_mc
 * 
 * \brief Namespace that contains mc package classes and variables.
 *
 * The rtt_mc namespace contains all classes and variables necessary to use
 * the mc package.  This namespace contains some "free" math functions and
 * constants in Constants.hh and Math.hh that are usefull in Monte Carlo
 * codes.

 * Additionally, this namespace holds all general purpose Monte Carlo
 * packages that exist in the mc package.  See \ref mc_overview for more
 * information.
 */
//===========================================================================//

#include <string>

namespace rtt_mc 
{
    const std::string release();
}

#endif                          // __mc_Release_hh__

//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
