//----------------------------------*-C++-*----------------------------------//
// P1Diffusion_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:00:28 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion.t.hh"
#include "SolverP1Diff.hh"

#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;

using namespace rtt_P1Diffusion;

typedef SolverP1Diff<MT> MS;

template
class P1Diffusion<MT,MS>;

//---------------------------------------------------------------------------//
//                              end of P1Diffusion_pt.cc
//---------------------------------------------------------------------------//
