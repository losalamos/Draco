//----------------------------------*-C++-*----------------------------------//
// Release.cc
// B.T. Adams
// Wed Jun 7 10:33:26 2000
/*! 
 * \file   meshReaders_Services/Release.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for meshReaders_Services library release 
 *         function.
 */
//---------------------------------------------------------------------------//
// @> Release function implementation for meshReaders_Servicess library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_meshReaders_Services
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form meshReaders_Servicess-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "meshReaders_Services(draco-4_3_0)";
    return pkg_release;
}

}  // end of rtt_meshReaders_Servicess namespace

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
