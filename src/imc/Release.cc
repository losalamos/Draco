//----------------------------------*-C++-*----------------------------------//
// Release.cc
// Thomas M. Evans
// Wed Apr 14 15:36:52 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for imc library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_imc
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form imc_#.#.# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)imc-1_3_0a";
    return pkg_release;
}

}  // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
