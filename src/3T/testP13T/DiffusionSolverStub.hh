//----------------------------------*-C++-*----------------------------------//
// DiffusionSolverStub.hh
// Randy M. Roberts
// Fri Mar 20 15:42:31 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_DiffusionSolverStub_hh__
#define __3T_testP13T_DiffusionSolverStub_hh__

#include "ds++/SP.hh"
#include <iostream>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class DiffusionSolverStub - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class DiffusionSolverStub
{
    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef typename MT::ccsf ccsf;
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::bsbf bsbf;
    typedef fcdsf FluxField;
    typedef fcdsf DiscFluxField;

    // DATA

  private:

    SP<MT> spmesh;
    
  public:

    // CREATORS

    DiffusionSolverStub(const SP<MT> &spmesh_) : spmesh(spmesh_) { }
    
    // MANIPULATORS
    
    // ACCESSORS

    const SP<MT> &getMesh() const { return spmesh; }

    void print(std::ostream &os) const
    {
	os << "in DiffusionSolverStub::print(), "
	   << "this: " << (void *)this
	   << " spmesh: ";
	spmesh->print(os);
    }
    
    void solve(const fcdsf &diffCoeff, const ccsf &removalCoeff,
	       const ccsf &source, const DiscFluxField &fluxSource,
	       const bsbf &boundary, ccsf &results,
	       FluxField &resultsFlux) const;

    void solve(const fcdsf &diffCoeff, const ccsf &removalCoeff,
	       const ccsf &source, const bsbf &boundary, ccsf &results) const;

  private:
    
    // IMPLEMENTATION
};

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_DiffusionSolverStub_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/DiffusionSolverStub.hh
//---------------------------------------------------------------------------//
