//----------------------------------*-C++-*----------------------------------//
// Version.cc
// Shawn Pautz
// Tue Mar 23 15:37:40 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Version function implementation for stopwatch library
//---------------------------------------------------------------------------//

#include "Version.hh"

namespace rtt_stopwatch
{
    using std::string;

  // function definition for Version, define the local version number for
  // this library in the form stopwatch_#.#.# in pkg_version variable
    const string version()
    {
	string pkg_version = "stopwatch_1.0.0";
	return pkg_version;
    }
}

//---------------------------------------------------------------------------//
//                              end of Version.cc
//---------------------------------------------------------------------------//
