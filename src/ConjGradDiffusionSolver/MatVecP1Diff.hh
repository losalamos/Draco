//----------------------------------*-C++-*----------------------------------//
// MatVecP1Diff.hh
// Randy M. Roberts
// Fri Sep 25 13:43:12 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_MatVecP1Diff_hh__
#define __ConjGradDiffusionSolver_MatVecP1Diff_hh__

#include "MatrixP1Diff.hh"
#include "ds++/SP.hh"

namespace rtt_ConjGradDiffusionSolver
{
 
// Since rtt_ConjGradDiffusionSolver is a limited namespace
// I risk putting in this using declaration for rtt_dsxx::SP,
// and rtt_dsxx::Mat1.
 
using rtt_dsxx::SP;

//===========================================================================//
// class MatVecP1Diff - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Matrix>
class MatVecP1Diff
{

    // NESTED CLASSES AND TYPEDEFS

  public:
     
    typedef typename Matrix::value_type value_type;
    typedef typename Matrix::ccsf ccsf;
     
  private:

    // DATA

    SP<const Matrix> spMatrix;
     
  public:

    // CREATORS
    
    MatVecP1Diff(const SP<const Matrix> &spMatrix_)
	: spMatrix(spMatrix_)
    {
	// empty
    }
     
    // MANIPULATORS
    
    // ACCESSORS

    ccsf operator()(const ccsf &x) const
    {
	ccsf b(x);
	spMatrix->multiply(b, x);
	return b;
    }

    // none

  private:
    
    // IMPLEMENTATION

    // none
};

} // namespace rtt_ConjGradDiffusionSolver

#endif                          // __ConjGradDiffusionSolver_MatVecP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGradDiffusionSolver/MatVecP1Diff.hh
//---------------------------------------------------------------------------//
