//----------------------------------*-C++-*----------------------------------//
// testP1Diffusion_pt.cc
// Randy M. Roberts
// Tue Sep 29 09:09:19 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testP1Diffusion.hh"
#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;

#include "P1Diffusion/SolverP1Diff.t.cc"
template class rtt_P1Diffusion::SolverP1Diff<MT>;

typedef rtt_P1Diffusion::SolverP1Diff<MT> Solver;

#include "P1Diffusion/P1Diffusion.t.cc"
template class rtt_P1Diffusion::P1Diffusion<MT, Solver >;

#include "P1Diffusion/MatVecP1Diff.t.cc"

template
class rtt_P1Diffusion::MatVecP1Diff<Solver::Matrix>;

#include "P1Diffusion/PreCondP1Diff.t.cc"

template
class rtt_P1Diffusion::PreCondP1Diff<Solver::Matrix>;

#include "P1Diffusion/MatrixP1Diff.t.cc"

template
class rtt_P1Diffusion::MatrixP1Diff<MT>;

//---------------------------------------------------------------------------//
//                              end of testP1Diffusion_pt.cc
//---------------------------------------------------------------------------//
