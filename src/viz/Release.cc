//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Release.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:29:46 2000
 * \brief  Release function implementation for viz library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_viz
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form viz-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "viz(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_viz

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
