//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/Release.cc
 * \author Randy M. Roberts
 * \date   Tue Jan 25 09:43:02 2000
 * \brief  Release function implementation for LAMG library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_LAMG
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form LAMG-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)LAMG-1_0_0";
    return pkg_release;
}

}  // end of rtt_LAMG

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
