//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Thomas M. Evans
// Tue Apr 13 17:58:54 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for mc library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_mc
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form mc_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "mc_1.0.0";
    return pkg_release;
}

}  // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
