//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   RTT_Format_Reader/Release.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_RTT_Format_Reader
{

using std::string;

// function definition for Release, define the local version number for
// this library in the form RTT_Format_Reader-#_#_# in pkg_version variable
const string release()
{
    string pkg_release = "RTT_Format_Reader(draco-3_0_0)";
    return pkg_release;
}

}  // end of rtt_RTT_Format_Reader namespace

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
