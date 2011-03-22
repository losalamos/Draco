//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Release.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  9 11:05:13 2005
 * \brief  Release function implementation for utils library
 * \note   Copyright (C) 2004-2006 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * 2010-11-29 This component was moved from clubimc/utils to
 * draco/diagnostics. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_diagnostics
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form diagnostics-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "diagnostics(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_diagnostics

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
