//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/PreComputedState_pt.cc
 * \author Randy M. Roberts
 * \date   Wed Jan 12 17:43:31 2000
 * \brief  Template instantiation file for PreComputedState.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../PreComputedState.t.hh"
#include "Mesh_XYZFactory.hh"

#include <iostream>

namespace rtt_LAMGDiffusionSolver
{

using rtt_LAMGDiffusionSolver_test::Mesh_XYZFactory;

typedef Mesh_XYZFactory::MT MT;

template
PreComputedState::PreComputedState<MT>(const MT::FieldConstructor &fCtor_in,
				       const MT &mesh);

} // end namespace rtt_LAMGDiffusionSolver_test

//---------------------------------------------------------------------------//
//                              end of PreComputedState_pt.cc
//---------------------------------------------------------------------------//
