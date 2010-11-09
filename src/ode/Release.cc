//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ode/Release.cc
 * \author Kent G. Budge
 * \date   Tue Nov  9 15:21:02 2010
 * \brief  Release function implementation for ode library
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_ode
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form ode-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "ode(draco-#_#_#)";
    return pkg_release;
}

}  // end of rtt_ode

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
