//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   P1Diffusion/test/DiffusionSelector.hh
 * \author Randy M. Roberts
 * \date   Tue Apr 25 08:14:28 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_test_DiffusionSelector_hh__
#define __P1Diffusion_test_DiffusionSelector_hh__

#include <P1Diffusion/config.h>

#ifdef SELECTOR_LAMG_TEST
#include "LAMGDiffusionSolver/SolverP1Diff.hh"
#endif

#ifdef SELECTOR_PCG_TEST
#include "PCGDiffusionSolver/SolverP1Diff.hh"
#include "PCGDiffusionSolver/pcg_DB.hh"
#endif

#ifdef SELECTOR_CONJGRAD_TEST
#include "ConjGradDiffusionSolver/SolverP1Diff.hh"
#endif

namespace rtt_P1Diffusion_test
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
#ifdef SELECTOR_LAMG_TEST
    typedef rtt_LAMGDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
#endif
  
#ifdef SELECTOR_PCG_TEST
    typedef rtt_PCGDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
    typedef rtt_PCGDiffusionSolver::pcg_DB Options;

    template<class NML_GROUP, class TDB>
    static Options create(const NML_GROUP &g, const TDB &tdb)
    {
	Options options("pcg");
	options.setup_namelist(g);
	return options;
    }
    
#endif
  
#ifdef SELECTOR_CONJGRAD_TEST
    typedef rtt_ConjGradDiffusionSolver::SolverP1Diff<MT> SolverP1Diff;
    typedef typename SolverP1Diff::Options Options;

    template<class NML_GROUP, class TDB>
    static Options create(const NML_GROUP &g, const TDB &tdb)
    {
	Options options(tdb.maxIterations, tdb.epsilon);
	return options;
    }
    
#endif
  
};

} // end namespace rtt_P1Diffusion_test

#endif                          // __P1Diffusion_test_DiffusionSelector_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/test/DiffusionSelector.hh
//---------------------------------------------------------------------------//
