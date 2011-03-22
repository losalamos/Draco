//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   traits/Release.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:56:25 2000
 * \brief  Release function implementation for traits library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_traits
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form traits-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "traits(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_traits

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
