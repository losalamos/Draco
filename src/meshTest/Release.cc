//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Release.cc
 * \author Randy M. Roberts
 * \date   Fri Aug 20 13:04:58 1999
 * \brief  Release function implementation for meshTest library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_meshTest
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form meshTest-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)meshTest-1_3_0";
    return pkg_release;
}

}  // end of rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
