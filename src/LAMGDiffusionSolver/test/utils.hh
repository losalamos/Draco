//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/utils.hh
 * \author Randy M. Roberts
 * \date   Wed Oct 27 14:32:09 1999
 * \brief  utilities for testing the LAMGDiffusionSolver
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_test_utils_hh__
#define __LAMGDiffusionSolver_test_utils_hh__

#include "Mesh_XYZFactory.hh"
#include <string>

namespace rtt_LAMGDiffusionSolver_test
{

void version(const std::string &progname);

Mesh_XYZFactory getMTFactory(const std::string &filename);

} // end namespace rtt_LAMGDiffusionSolver_test

#endif                          // __LAMGDiffusionSolver_test_utils_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/test/utils.hh
//---------------------------------------------------------------------------//
