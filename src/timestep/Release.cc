//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   timestep/Release.cc
 * \author Thomas M. Evans
 * \date   Mon Apr 19 21:36:00 2004
 * \brief  Release function implementation for timestep library
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_timestep
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form timestep-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "timestep(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_timestep

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
