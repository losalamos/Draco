//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   P1Diffusion/DiffusionSelector.hh
 * \author Randy M. Roberts
 * \date   Tue Apr 25 08:14:28 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_DiffusionSelector_hh__
#define __P1Diffusion_DiffusionSelector_hh__

#include <P1Diffusion/config.h>

#ifdef SELECTOR_LAMG
#include "LAMGDiffusionSolver/SolverP1Diff.hh"
#endif

#ifdef SELECTOR_PCG
#include "PCGDiffusionSolver/SolverP1Diff.hh"
#include "PCGDiffusionSolver/pcg_DB.hh"
#endif

#ifdef SELECTOR_CONJGRAD
#include "ConjGradDiffusionSolver/SolverP1Diff.hh"
#endif

namespace rtt_P1Diffusion
{
 
//===========================================================================//
/*!
 * \class DiffusionSelector
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
struct DiffusionSelector 
{
#ifdef SELECTOR_LAMG
    typedef rtt_LAMGDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
#endif
  
#ifdef SELECTOR_PCG
    typedef rtt_PCGDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
    typedef rtt_PCGDiffusionSolver::pcg_DB Options;
#endif
  
#ifdef SELECTOR_CONJGRAD
    typedef rtt_ConjGradDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
    typedef typename SolverP1Diff::Options Options;
#endif
  
};

} // end namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_DiffusionSelector_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/DiffusionSelector.hh
//---------------------------------------------------------------------------//
