//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/Release.cc
 * \author Randy M. Roberts
 * \date   Mon Apr 24 08:14:47 2000
 * \brief  Release function implementation for ConjGrad library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_ConjGrad
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form ConjGrad-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)ConjGrad-1_0_0";
    return pkg_release;
}

}  // end of rtt_ConjGrad

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
