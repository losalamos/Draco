//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for POOMA_MT test suite
//---------------------------------------------------------------------------//

#include "TestRelease.hh"

namespace rtt_POOMA_MT_test
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form POOMA_MT_test-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)POOMA_MT_test-1_1_0";
    return pkg_release;
}

}  // end of rtt_POOMA_MT_test

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
