//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/Release.cc
 * \author Randy M. Roberts
 * \date   Mon Apr 24 15:22:41 2000
 * \brief  Release function implementation for ConjGradDiffusionSolver library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_ConjGradDiffusionSolver
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form ConjGradDiffusionSolver-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)ConjGradDiffusionSolver-1_0_0";
    return pkg_release;
}

}  // end of rtt_ConjGradDiffusionSolver

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
