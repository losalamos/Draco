//----------------------------------*-C++-*----------------------------------//
// Release.cc
// John Gulick
// Tue Aug 24 13:08:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function implementation for fourier library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_fourier
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form fourier-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)fourier-1_0_0";
    return pkg_release;
}

}  // end of rtt_fourier

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
