//----------------------------------*-C++-*----------------------------------//
// DiffusionSolverStub.cc
// Randy M. Roberts
// Fri Mar 20 15:42:31 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/DiffusionSolverStub.hh"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
template<class MT>
void DiffusionSolverStub<MT>::solve(const fcdsf &diffCoeff,
				    const ccsf &removalCoeff,
				    const ccsf &source,
				    const DiscFluxField &fluxSource,
				    const bsbf &boundary,
				    ccsf &results,
				    FluxField &resultsFlux) const
{
    results = source / removalCoeff;
    resultsFlux = 0.0;
}

template<class MT>
void DiffusionSolverStub<MT>::solve(const fcdsf &diffCoeff,
				    const ccsf &removalCoeff,
				    const ccsf &source,
				    ccsf &results) const
{
    results = source / removalCoeff;
}


END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of DiffusionSolverStub.cc
//---------------------------------------------------------------------------//
