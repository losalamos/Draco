//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Release.cc
 * \author Thomas M. Evans
 * \date   Thu Nov 18 11:48:00 1999
 * \brief  Release function implementation for mc library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_mc
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form mc-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)mc-1_1_0x";
    return pkg_release;
}

}  // end of rtt_mc

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
