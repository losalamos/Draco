//----------------------------------*-C++-*----------------------------------//
// Release.cc
// <user>
// <date>
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for <pkg> library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_<spkg>
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form <spkg>_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)<spkg>_<start>#.#.#";
    return pkg_release;
}

}  // end of rtt_<spkg>

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
