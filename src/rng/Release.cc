//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Thomas M. Evans
// Thu May 27 15:24:02 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for rng library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_rng
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form rng_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)rng-1_0_0";
    return pkg_release;
}

}  // end of rtt_rng

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
