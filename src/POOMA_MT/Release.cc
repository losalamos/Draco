//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Sun Aug 22 08:26:12 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for POOMA_MT library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_POOMA_MT
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form POOMA_MT-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)POOMA_MT-1_1_0";
    return pkg_release;
}

}  // end of rtt_POOMA_MT

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
