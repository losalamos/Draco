//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/Release.cc
 * \author Randy M. Roberts
 * \date   Wed Feb  9 10:57:46 2000
 * \brief  Release function implementation for units library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_units
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form units-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "units(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_units

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
