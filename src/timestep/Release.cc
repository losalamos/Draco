//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Fri Aug 27 10:33:26 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for timestep library
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_timestep
{

// function definition for Release, define the local version number for
// this library in the form timestep-#_#_# in pkg_version variable
const std::string release()
{
    std::string pkg_release = "@(#)timestep-1_0_2";
    return pkg_release;
}

}  // end of rtt_timestep

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
