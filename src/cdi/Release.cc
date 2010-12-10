//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/Release.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:18:14 2000
 * \brief  Release function implementation for cdi library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_cdi
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form cdi-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "cdi(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_cdi

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
