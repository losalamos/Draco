//----------------------------------*-C++-*----------------------------------//
// MatrixP1DiffTraits.hh
// Randy M. Roberts
// Tue Jun  8 09:50:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_MatrixP1DiffTraits_hh__
#define __LAMGDiffusionSolver_MatrixP1DiffTraits_hh__

#include "traits/MatrixFactoryTraits.hh"
#include "MatrixP1Diff.hh"
#include "PreComputedState.hh"
#include "LAMG/CompressedRowStorage.hh"
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

template<>
class MatrixFactoryTraits<rtt_LAMGDiffusionSolver::MatrixP1Diff>
{

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_LAMGDiffusionSolver::MatrixP1Diff Matrix;

  public:

    typedef rtt_LAMGDiffusionSolver::PreComputedState PreComputedState;
    
  private:

    // DATA
    
  public:

    // STATIC CLASS METHODS

    template<class MT>
    static
    PreComputedState preComputeState(const typename MT::FieldConstructor &fCtor,
				     const MT &mesh)
    {
	// Return the needed precomputed state.
	
	return PreComputedState(fCtor, mesh);
    }

    template<class MT>
    static Matrix *create(const rtt_diffusion::P1Matrix<MT> &rep,
			  const PreComputedState &state)
    {
	typedef rtt_LAMG::CompressedRowStorage CRS;
	CRS crs(state.rowPointer(), state.colIndex(),
		state.val(rep.diagonal().begin(), rep.diagonal().end(),
			  rep.offDiagonal().begin(), rep.offDiagonal().end()));
	return new Matrix(crs);
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

#endif                          // __LAMGDiffusionSolver_MatrixP1DiffTraits_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/MatrixP1DiffTraits.hh
//---------------------------------------------------------------------------//
