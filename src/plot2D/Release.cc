//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   plot2D/Release.cc
 * \author Rob Lowrie
 * \date   Mon Apr 15 10:05:15 2002
 * \brief  Release function implementation for plot2D library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_plot2D
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form plot2D-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "plot2D(draco-4_3_0)";
    return pkg_release;
}

}  // end of rtt_plot2D

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
