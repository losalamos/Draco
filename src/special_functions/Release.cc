//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/Release.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 09:48:39 2000
 * \brief  Release function implementation for sf library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_sf
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form sf-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "sf(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_sf

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
