//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   shared_lib/Release.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 21 14:30:27 2004
 * \brief  Release function implementation for shared_lib library
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_shared_lib
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form shared_lib-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "shared_lib(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_shared_lib

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
