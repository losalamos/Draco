//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/TestMatrixGen_pt.cc
 * \author Randy M. Roberts
 * \date   Wed Jan 12 17:27:35 2000
 * \brief  Template instantiation file for TestMatrixGen.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMatrixGen.t.hh"
#include "Mesh_XYZFactory.hh"

namespace rtt_LAMGDiffusionSolver_test
{

template class TestMatrixGen<Mesh_XYZFactory>;

} // end namespace rtt_LAMGDiffusionSolver_test

//---------------------------------------------------------------------------//
//                              end of TestMatrixGen_pt.cc
//---------------------------------------------------------------------------//
