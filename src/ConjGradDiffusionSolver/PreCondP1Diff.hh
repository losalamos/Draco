//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/PreCondP1Diff.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 16:10:39 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_PreCondP1Diff_hh__
#define __ConjGradDiffusionSolver_PreCondP1Diff_hh__

#include "MatrixP1Diff.hh"
#include "ds++/SP.hh"

namespace rtt_ConjGradDiffusionSolver
{

// Since rtt_ConjGradDiffusionSolver is a limited namespace
// I risk putting in this using declaration for rtt_dsxx::SP,
// and rtt_dsxx::Mat1.
 
using rtt_dsxx::SP;
 
//===========================================================================//
/*!
 * \class PreCondP1Diff
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Matrix>
class PreCondP1Diff 
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
    
    PreCondP1Diff(const SP<const Matrix> &spMatrix_)
	: spMatrix(spMatrix_)
    {
	// empty
    }

    //Defaulted: PreCondP1Diff(const PreCondP1Diff &rhs);
    //Defaulted: ~PreCondP1Diff();

    // MANIPULATORS
    
    //Defaulted: PreCondP1Diff& operator=(const PreCondP1Diff &rhs);

    // ACCESSORS

    ccsf operator()(const ccsf &x) const
    {
	ccsf b(x);
	spMatrix->jacobiIteration(b, x);
	return b;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGradDiffusionSolver

#endif                          // __ConjGradDiffusionSolver_PreCondP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGradDiffusionSolver/PreCondP1Diff.hh
//---------------------------------------------------------------------------//
