//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for mesh library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_mesh
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form mesh-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)mesh-1_8_0";
    return pkg_release;
}

}  // end of rtt_mesh

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
