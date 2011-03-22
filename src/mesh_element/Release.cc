//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_element/Release.cc
 * \author Kelly Thompson
 * \date   Mon May 24 16:54:16 2004
 * \brief  Release function implementation for mesh_element library
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_mesh_element
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form mesh_element-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "mesh_element(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_mesh_element

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
