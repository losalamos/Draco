//----------------------------------*-C++-*----------------------------------//
// Release.cc
// B.T. Adams
// Wed Jun 7 10:33:26 2000
/*! 
 * \file   RTT_Format_Reader/Release.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Readers library release function.
 */
//---------------------------------------------------------------------------//
// @> Release function implementation for RTT_Format_Readers library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_RTT_Format_Reader
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form RTT_Format_Readers-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "@(#)RTT_Format_Reader-1_0_0";
    return pkg_release;
}

}  // end of rtt_RTT_Format_Readers namespace

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
