//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Fri Aug 20 13:04:58 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for meshTest library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_meshTest
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form meshTest-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)meshTest-1_0_0";
    return pkg_release;
}

}  // end of rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
