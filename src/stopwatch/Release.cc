//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Shawn Pautz
// Tue Mar 23 15:37:40 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for stopwatch library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_stopwatch
{
    using std::string;

  // function definition for Release, define the local version number for
  // this library in the form stopwatch_#.#.# in pkg_version variable
    const string release()
    {
	string pkg_release = "stopwatch_1.0.0";
	return pkg_release;
    }
}  // end of rtt_stopwatch

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
