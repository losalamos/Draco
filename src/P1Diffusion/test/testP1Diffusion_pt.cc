//----------------------------------*-C++-*----------------------------------//
// testP1Diffusion_pt.cc
// Randy M. Roberts
// Tue Sep 29 09:09:19 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testP1Diffusion.hh"
#include "DiffusionSelector.hh"
#include "mesh/Mesh_XYZ.hh"

typedef Mesh_XYZ MT;
typedef rtt_P1Diffusion_test::DiffusionSelector<MT>::SolverP1Diff Solver;

#ifdef SELECTOR_PCG_TEST

#include "PCGDiffusionSolver/SolverP1Diff.t.hh"

template class rtt_PCGDiffusionSolver::SolverP1Diff<MT>;

#include "../P1Diffusion.t.hh"
#include "PCGDiffusionSolver/MatrixP1DiffTraits.hh"
template class rtt_P1Diffusion::P1Diffusion<MT, Solver>;

#include "PCGDiffusionSolver/MatVecP1Diff.t.hh"

template
class rtt_PCGDiffusionSolver::MatVecP1Diff<Solver::Matrix>;

#include "PCGDiffusionSolver/PreCondP1Diff.t.hh"

template
class rtt_PCGDiffusionSolver::PreCondP1Diff<Solver::Matrix>;

#include "PCGDiffusionSolver/MatrixP1Diff.t.hh"

template
class rtt_PCGDiffusionSolver::MatrixP1Diff<MT>;

#endif

#ifdef SELECTOR_CONJGRAD_TEST

#include "../P1Diffusion.t.hh"
template class rtt_P1Diffusion::P1Diffusion<MT, Solver>;

#endif

//---------------------------------------------------------------------------//
//                              end of testP1Diffusion_pt.cc
//---------------------------------------------------------------------------//
