//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/Release.cc
 * \author Randy M. Roberts
 * \date   Mon Jan 10 15:40:00 2000
 * \brief  Release function implementation for LAMGDiffusionSolver library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_LAMGDiffusionSolver
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form LAMGDiffusionSolver-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)LAMGDiffusionSolver-1_1_0";
    return pkg_release;
}

}  // end of rtt_LAMGDiffusionSolver

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
