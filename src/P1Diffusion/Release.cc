//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   P1Diffusion/Release.cc
 * \author Randy M. Roberts
 * \date   Mon Jan  3 12:50:18 2000
 * \brief  Release function implementation for P1Diffusion library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_P1Diffusion
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form P1Diffusion-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)P1Diffusion-1_0_0";
    return pkg_release;
}

}  // end of rtt_P1Diffusion

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
