//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders/Release.cc
 * \author B.T. Adams
 * \date   Wed Apr 14 10:33:26 1999
 * \brief  Implementation file for meshReaders library release function.
 */
//---------------------------------------------------------------------------//
// @> Release function implementation for meshReaders library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_meshReaders
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form meshReader-#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "meshReaders(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_meshReaders namespace

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
