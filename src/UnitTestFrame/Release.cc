//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/Release.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 10:10:23 2000
 * \brief  Release function implementation for UnitTestFrame library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_UnitTestFrame
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form UnitTestFrame-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)UnitTestFrame-1_0_0";
    return pkg_release;
}

}  // end of rtt_UnitTestFrame

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
