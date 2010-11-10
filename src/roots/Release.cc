//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   roots/Release.cc
 * \author Kent G. Budge
 * \date   Tue Nov  9 16:49:40 2010
 * \brief  Release function implementation for roots library
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_roots
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form roots-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "roots(draco-#_#_#)";
    return pkg_release;
}

}  // end of rtt_roots

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
