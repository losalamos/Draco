//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   radphys/Release.cc
 * \author Randy M. Roberts
 * \date   Tue Feb  1 09:30:08 2000
 * \brief  Release function implementation for radphys library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_radphys
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form radphys-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)radphys-1_3_0";
    return pkg_release;
}

}  // end of rtt_radphys

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
