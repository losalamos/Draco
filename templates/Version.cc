//----------------------------------*-C++-*----------------------------------//
// Version.cc
// <user>
// <date>
// $Id$
//---------------------------------------------------------------------------//
// @> Version function implementation for <pkg> library
//---------------------------------------------------------------------------//

#include "Version.hh"

namespace rtt_<spkg>
{
  // function definition for Version, define the local version number for
  // this library in the form <spkg>_#.#.# in pkg_version variable
    const string version()
    {
	string pkg_version = "<spkg>_<start>#.#.#";
	return pkg_version;
    }
}

//---------------------------------------------------------------------------//
//                              end of Version.cc
//---------------------------------------------------------------------------//
