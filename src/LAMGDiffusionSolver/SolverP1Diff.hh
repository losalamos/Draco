//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/SolverP1Diff.hh
 * \author Randy M. Roberts
 * \date   Tue Jan 25 08:40:08 2000
 * \brief  This is the LAMG-based MatrixSolver to be used within P1Diffusion.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_SolverP1Diff_hh__
#define __LAMGDiffusionSolver_SolverP1Diff_hh__

#include "LAMG/Solver.hh"
#include "ds++/SP.hh"

namespace rtt_LAMGDiffusionSolver
{

// Forward Declarations

class MatrixP1Diff;

// Inside a small namespace I feel I can put in a using declaration.
 
 using dsxx::SP;

//===========================================================================//
/*!
 * \class SolverP1Diff
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class SolverP1Diff 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    rtt_LAMG::Solver linsolver;
    
  public:

    // CREATORS
    
    SolverP1Diff(const rtt_LAMG::Options& opts)
    : linsolver(opts)
    {
	// empty
    }
    
    //DEFAULTED: SolverP1Diff(const SolverP1Diff &rhs);
    //DEFAULTED: ~SolverP1Diff();

    // MANIPULATORS

    template<class FT>
    void solve(FT &x, const SP<const MatrixP1Diff> &spMatrix,
	       const FT &b)
    {
	Require(x.size() == b.size());
	linsolver.solve(x, spMatrix->lamgMatrixDcsrR(), b);
    }

    //DEFAULTED: SolverP1Diff& operator=(const SolverP1Diff &rhs);

    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_LAMGDiffusionSolver

#endif                          // __LAMGDiffusionSolver_SolverP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/SolverP1Diff.hh
//---------------------------------------------------------------------------//
