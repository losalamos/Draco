//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diffusion/Release.cc
 * \author Randy M. Roberts
 * \date   Tue Feb  1 09:24:50 2000
 * \brief  Release function implementation for diffusion library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_diffusion
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form diffusion-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)diffusion-1_3_0";
    return pkg_release;
}

}  // end of rtt_diffusion

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
