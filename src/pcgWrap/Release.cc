//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   pcgWrap/Release.cc
 * \author Rob Lowrie
 * \date   Thu Jul 15 09:44:12 1999
 * \brief  Release function implementation for pcgWrap library
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_pcgWrap
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form pcgWrap-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "pcgWrap(draco-5_2_0)";
    return pkg_release;
}

}  // end of rtt_pcgWrap

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
