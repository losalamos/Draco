//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/Release.cc
 * \author <user>
 * \date   <date>
 * \brief  Release function implementation for <pkg> library
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace <namespace>
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form <spkg>-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "<spkg>(draco-<start>#_#_#)";
    return pkg_release;
}

}  // end of <namespace>

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
