//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PCGDiffusionSolver/Release.cc
 * \author Randy M. Roberts
 * \date   Mon Jan  3 12:46:57 2000
 * \brief  Release function implementation for PCGDiffusionSolver library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_PCGDiffusionSolver
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form PCGDiffusionSolver-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)PCGDiffusionSolver-1_2_0";
    return pkg_release;
}

}  // end of rtt_PCGDiffusionSolver

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
