//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/Release.cc
 * \author Thomas M. Evans
 * \date   Mon Apr 19 22:08:44 2004
 * \brief  Release function implementation for RTT_Format_Reader library
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_RTT_Format_Reader
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form RTT_Format_Reader-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "RTT_Format_Reader(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
