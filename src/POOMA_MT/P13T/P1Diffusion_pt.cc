//----------------------------------*-C++-*----------------------------------//
// P1Diffusion_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:00:28 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/P1Diffusion.t.cc"
#include "P1Diffusion/SolverP1Diff.hh"

#include "POOMA_MT/PoomaMesh_XYZ.hh"

typedef UniformCartesian<3> PoomaMesh_t;
typedef PoomaMesh_XYZ<PoomaMesh_t> MT;

using namespace rtt_P1Diffusion;

typedef SolverP1Diff<MT> MS;

template
class P1Diffusion<MT,MS>;

//---------------------------------------------------------------------------//
//                              end of P1Diffusion_pt.cc
//---------------------------------------------------------------------------//
