//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Thu Jul 15 09:41:09 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for c4 library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_c4
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form c4-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)c4-1_0_0";
    return pkg_release;
}

}  // end of rtt_c4

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
