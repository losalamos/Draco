//----------------------------------*-C++-*----------------------------------//
// SolverP1Diff_pt.cc
// Randy M. Roberts
// Thu Oct  1 16:02:32 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/SolverP1Diff.t.cc"
#include "POOMA_MT/PoomaMesh_XYZ.hh"

typedef UniformCartesian<3> PoomaMesh_t;
typedef PoomaMesh_XYZ<PoomaMesh_t> MT;

using namespace rtt_P1Diffusion;

template
class SolverP1Diff<MT>;

//---------------------------------------------------------------------------//
//                              end of SolverP1Diff_pt.cc
//---------------------------------------------------------------------------//
