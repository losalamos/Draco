//----------------------------------*-C++-*----------------------------------//
// Release.cc
// B.T. Adams
// Wed Apr 14 15:36:52 1999
// $Id$
/*! 
 * \file   amr_mesh/Release.cc
 * \author B.T. Adams
 * \date   Wed Apr 14 10:33:26 1999
 * \brief  Implementation file for amr_mesh library release function.
 */
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
