//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nml/Release.cc
 * \author Randy M. Roberts
 * \date   Wed Feb  9 10:52:51 2000
 * \brief  Release function implementation for nml library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_nml
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form nml-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)nml-1_1_0";
    return pkg_release;
}

}  // end of rtt_nml

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
