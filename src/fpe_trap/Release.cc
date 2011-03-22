//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/Release.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:00:50 2005
 * \brief  Release function implementation for fpe_trap library
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_fpe_trap
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form fpe_trap-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "fpe_trap(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_fpe_trap

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
