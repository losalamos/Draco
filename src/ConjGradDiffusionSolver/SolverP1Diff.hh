//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/SolverP1Diff.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 15:32:40 2000
 * \brief  This is the ConjGrad-based MatrixSolver to be used within P1Diffusion.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_SolverP1Diff_hh__
#define __ConjGradDiffusionSolver_SolverP1Diff_hh__

#include "MatrixP1Diff.hh"
#include "MatVecP1Diff.hh"
#include "PreCondP1Diff.hh"

#include "ConjGrad/ConjGrad.hh"
#include "ConjGrad/ConjGradTraits.hh"
#include "ds++/SP.hh"

#include <iostream>

namespace rtt_ConjGradDiffusionSolver
{
 
// Inside a small namespace I feel I can put in a using declaration.
 
 using rtt_dsxx::SP;

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

template<class MT>
class SolverP1Diff 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef typename MT::ccsf ccsf;
    typedef typename MT::FieldConstructor FieldConstructor;
    typedef MatrixP1Diff<MT> Matrix;
    typedef MatVecP1Diff<Matrix> MatVec;
    typedef PreCondP1Diff<Matrix> PreCond;
    typedef typename ccsf::value_type value_type;

    struct Options
    {
	int maxIters;
	double eps;
        bool usePreconditioner;
	bool verbose;
	Options(int maxIters_in, double eps_in, bool usePreconditioner_in,
                bool verbose_in)
	    : maxIters(maxIters_in), eps(eps_in),
              usePreconditioner(usePreconditioner_in), verbose(verbose_in)
	{
	    // empty
	}
    };

  private:

    // DATA
    
    SP<const MT> spMesh;
    FieldConstructor fCtor;
    Options options;
    
  public:

    // CREATORS
    
    SolverP1Diff(const SP<const MT>& spMesh_in,
		 const FieldConstructor &fCtor_in,
		 int maxIters_in, double eps_in,
                 bool usePreconditioner_in,
                 bool verbose_in)
	: spMesh(spMesh_in), fCtor(fCtor_in), options(maxIters_in, eps_in,
                                                      usePreconditioner_in,
						      verbose_in)
    {
	// empty
    }
    
    SolverP1Diff(const SP<const MT>& spMesh_in,
		 const FieldConstructor &fCtor_in,
		 const Options &options_in)
	: spMesh(spMesh_in), fCtor(fCtor_in), options(options_in)
    {
	// empty
    }
    
    SolverP1Diff(const SP<const MT>& spMesh_in,
		 const Options &options_in)
	: spMesh(spMesh_in), fCtor(spMesh_in), options(options_in)
    {
	// empty
    }
    
    //Defaulted: SolverP1Diff(const SolverP1Diff &rhs);
    //Defaulted: ~SolverP1Diff();

    // MANIPULATORS
    
    //Defaulted: SolverP1Diff& operator=(const SolverP1Diff &rhs);

    void solve(ccsf &phi, const SP<const Matrix> &spMatrix,
	       const ccsf &brhs)
    {
	int iter;
	ccsf r(fCtor);
	using rtt_ConjGrad::conjGrad;

        if (options.usePreconditioner)
            conjGrad(phi, iter, brhs, MatVec(spMatrix), options.maxIters,
                     options.eps, PreCond(spMatrix), r);
        else
            conjGrad(phi, iter, brhs, MatVec(spMatrix), options.maxIters,
                     options.eps, r);
        
	if (options.verbose)
	    std::cout << "SolverP1Diff: " << iter << " iterations, "
		      << "||r||: "
		      << rtt_ConjGrad::ConjGradTraits<ccsf>::Norm()(r)
		      << ", ||b||: "
		      << rtt_ConjGrad::ConjGradTraits<ccsf>::Norm()(brhs)
		      << std::endl;
    }

    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGradDiffusionSolver

#endif                          // __ConjGradDiffusionSolver_SolverP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGradDiffusionSolver/SolverP1Diff.hh
//---------------------------------------------------------------------------//
