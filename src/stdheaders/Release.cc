//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   stdheaders/Release.cc
 * \author Michelle L. Murillo
 * \date   Wed Jun 21 17:46:55 2000
 * \brief  Release function implementation for stdheaders library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_stdheaders
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form stdheaders-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "stdheaders(draco-4_0_0)";
    return pkg_release;
}

}  // end of rtt_stdheaders

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
