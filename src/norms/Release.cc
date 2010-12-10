//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/Release.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 09:06:57 2005
 * \brief  Release function implementation for norms library
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_norms
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form norms-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "norms(draco-6_0_0)";
    return pkg_release;
}

}  // end of rtt_norms

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
