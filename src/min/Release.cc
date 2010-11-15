//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   min/Release.cc
 * \author Kent G. Budge
 * \date   Mon Nov 15 10:20:40 2010
 * \brief  Release function implementation for min library
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_min
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form min-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "min(draco-#_#_#)";
    return pkg_release;
}

}  // end of rtt_min

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
