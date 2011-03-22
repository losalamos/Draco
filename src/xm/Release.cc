//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xm/Release.cc
 * \author Randy M. Roberts
 * \date   Wed Feb  9 10:59:38 2000
 * \brief  Release function implementation for xm library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_xm
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form xm-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "xm(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_xm

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
