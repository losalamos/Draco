//----------------------------------*-C++-*----------------------------------//
/*!
 * \file ds++/Release.cc
 * \author Randy M. Roberts
 * \date Thu Jul 15 09:31:44 1999
 * \brief Provides the function definition for Release.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_dsxx
{

// function definition for Release, define the local version number for
// this library in the form ds_#.#.# in pkg_version variable
const std::string release()
{
    std::string pkg_release = "@(#)ds-1_0_0";
    return pkg_release;
}

}  // end of rtt_dsxx

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
