//----------------------------------*-C++-*----------------------------------//
// MatrixP1DiffTraits.hh
// Randy M. Roberts
// Tue Jun  8 09:50:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __PCGDiffusionSolver_MatrixP1DiffTraits_hh__
#define __PCGDiffusionSolver_MatrixP1DiffTraits_hh__

#include "traits/MatrixFactoryTraits.hh"
#include "PCGDiffusionSolver/MatrixP1Diff.hh"
#include "diffusion/P1Matrix.hh"

namespace rtt_traits
{
 
//===========================================================================//
// class MatrixP1DiffTraits - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class MatrixFactoryTraits<rtt_PCGDiffusionSolver::MatrixP1Diff<MT> >
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::FieldConstructor FieldConstructor;
    
    typedef rtt_PCGDiffusionSolver::MatrixP1Diff<MT> Matrix;

  public:

    typedef int PreComputedState;
    
  private:

    // DATA
    
  public:

    // STATIC CLASS METHODS

    static PreComputedState preComputeState(const FieldConstructor &fCtor,
					    const MT &mesh)
    {
	// There is no needed precomputed state.
	
	return 42;
    }
    
    static Matrix *create(const rtt_diffusion::P1Matrix<MT> &rep,
			  const PreComputedState &state)
    {
	return new Matrix(rep.fCtor(), rep.diagonal(), rep.offDiagonal());
    }

    // CREATORS
    
    // ** none **

    // MANIPULATORS
    
    // ** none **

    // ACCESSORS

  private:
    
    // IMPLEMENTATION

    // ** none **
};

} // end namespace rtt_traits

#endif                          // __PCGDiffusionSolver_MatrixP1DiffTraits_hh__

//---------------------------------------------------------------------------//
//                              end of PCGDiffusionSolver/MatrixP1DiffTraits.hh
//---------------------------------------------------------------------------//
