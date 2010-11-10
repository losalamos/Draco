//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/Release.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 08:02:06 2010
 * \brief  Release function implementation for linear library
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_linear
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form linear-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "linear(draco-#_#_#)";
    return pkg_release;
}

}  // end of rtt_linear

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
