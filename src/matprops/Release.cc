//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   matprops/Release.cc
 * \author Randy M. Roberts
 * \date   Tue Nov 16 08:51:10 1999
 * \brief  Release function implementation for matprops library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_matprops
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form matprops-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)matprops-1_2_0";
    return pkg_release;
}

}  // end of rtt_matprops

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
