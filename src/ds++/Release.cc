//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Randy M. Roberts
// Thu Jul 15 09:31:44 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for ds++ library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_ds
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form ds_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)ds-1_0_0";
    return pkg_release;
}

}  // end of rtt_ds

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
