//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff.hh
// Randy M. Roberts
// Mon Sep 28 08:33:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_MatrixP1Diff_hh__
#define __ConjGradDiffusionSolver_MatrixP1Diff_hh__

#include "ds++/SP.hh"

namespace rtt_ConjGradDiffusionSolver
{

 // Since rtt_ConjGradDiffusionSolver is a limited namespace
 // I risk putting in this using declaration for rtt_dsxx::SP.
 
 using rtt_dsxx::SP;
 
 //===========================================================================//
 // class MatrixP1Diff - 
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
 class MatrixP1Diff
 {
     // NESTED CLASSES AND TYPEDEFS

   public:

     typedef typename MT::ccsf ccsf;
     typedef typename MT::fcdsf fcdsf;
     typedef typename MT::FieldConstructor FieldConstructor;

     typedef typename ccsf::value_type value_type;
     
     // DATA

   private:
     
     FieldConstructor fCtor;
     ccsf diagonal_m;
     fcdsf offDiagonal_m;

     // Mutable temporary data

     mutable fcdsf xFC;
     mutable fcdsf xSwap;
     mutable fcdsf bFC;
     mutable ccsf btmp;
     mutable ccsf xtmp;
     
   public:

     // CREATORS
    
     MatrixP1Diff(const FieldConstructor &fCtor_, const ccsf &diagonal_,
		  const fcdsf &offDiagonal_)
	 : fCtor(fCtor_), diagonal_m(diagonal_), offDiagonal_m(offDiagonal_),
	   xFC(fCtor_), xSwap(fCtor_), bFC(fCtor_), btmp(fCtor_), xtmp(fCtor_)
     {
	 // empty
     }
     
     // MANIPULATORS

     // ACCESSORS

     void multiply(ccsf &b, const ccsf &x) const
     {
	 b = diagonal_m * x;
	 
	 MT::gather(xFC, x, MT::OpAssign());
	 MT::swap_faces(xSwap, xFC, 0.0);

	 bFC = offDiagonal_m * xSwap;
	 
	 MT::scatter(b, bFC, MT::OpAddAssign());
     }

     void jacobiIteration(ccsf &b, const ccsf &x) const
     {
	 b = x / diagonal_m;
     }

     const ccsf &diagonal() const
     {
	 return diagonal_m;
     }

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_ConjGradDiffusionSolver

// Let's always include the traits class with us.

#include "MatrixP1DiffTraits.hh"

#endif                          // __ConjGradDiffusionSolver_MatrixP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGradDiffusionSolver/MatrixP1Diff.hh
//---------------------------------------------------------------------------//
