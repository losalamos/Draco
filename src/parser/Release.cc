//----------------------------------*-C++-*----------------------------------//
/*!
 * \file parser/Release.cc
 * \author Kent G. Budge
 * \date Thu Jul 15 09:31:44 1999
 * \brief Provides the function definition for Release.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_parser
{

// function definition for Release, define the local version number for
// this library in the form parser_#.#.# in pkg_version variable
const std::string release()
{
    std::string pkg_release = "parser(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_parser

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
