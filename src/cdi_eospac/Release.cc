//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/Release.cc
 * \author Kelly Thompson
 * \date   Fri Mar 30 10:43:31 2001
 * \brief  Release function implementation for cdi_eospac library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_cdi_eospac
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form cdi_eospac-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "cdi_eospac(draco-3_0_0)";
    return pkg_release;
}

}  // end of rtt_cdi_eospac

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
