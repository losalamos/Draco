//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   amr_mesh/Release.cc
 * \author Todd Adams
 * \date   Tue Sep 21 15:10:00 1999
 * \brief  Release function implementation for amr_mesh library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_amr_mesh
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form amr_mesh-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)amr_mesh-0_0_0";
    return pkg_release;
}

}  // end of rtt_amr_mesh

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
