//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/Release.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:18:14 2000
 * \brief  Release function implementation for cdi_gandolf library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_cdi_gandolf
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
    string pkg_release = "@(#)cdi-#_#_#";
    return pkg_release;
}

}  // end of namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
