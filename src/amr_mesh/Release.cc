//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Thomas M. Evans
// Wed Apr 14 15:36:52 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for amr library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_amr
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form amr_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)amr-1_0_0";
    return pkg_release;
}

}  // end of rtt_amr

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
