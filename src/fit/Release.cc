//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fit/Release.cc
 * \author Kent G. Budge
 * \date   Mon Nov 15 10:27:30 2010
 * \brief  Release function implementation for fit library
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_fit
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form fit-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "fit(draco-#_#_#)";
    return pkg_release;
}

}  // end of rtt_fit

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
