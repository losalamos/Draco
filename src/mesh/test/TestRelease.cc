//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for mesh test suite
//---------------------------------------------------------------------------//

#include "TestRelease.hh"

namespace rtt_mesh_test
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form mesh_test-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)mesh_test-1_2_0";
    return pkg_release;
}

}  // end of rtt_mesh_test

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
