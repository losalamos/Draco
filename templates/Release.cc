//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/Release.cc
 * \author <user>
 * \date   <date>
 * \brief  Release function implementation for <pkg> library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_<spkg>
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form <spkg>-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "<spkg>(draco-<start>#_#_#)";
    return pkg_release;
}

}  // end of rtt_<spkg>

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
