//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff.hh
// Randy M. Roberts
// Mon Sep 28 08:33:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __PCGDiffusionSolver_MatrixP1Diff_hh__
#define __PCGDiffusionSolver_MatrixP1Diff_hh__

#include "ds++/SP.hh"

namespace rtt_PCGDiffusionSolver
{

 // Since rtt_PCGDiffusionSolver is a limited namespace
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
     typedef fcdsf Vector;
     
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

     void multiply(ccsf &b, const ccsf &x) const;

     void multiply(const SP<ccsf> &b, const SP<const ccsf> &x) const
     {
	 multiply(*b, *x);
     }

     template<class FT>
     void multiply(FT &b, const FT &x) const;

     template<class FT>
     void multiply(const SP<FT> &b, const SP<const FT> &x) const
     {
	 multiply(*b, *x);
     }

     void jacobiIteration(ccsf &b, const ccsf &x) const;

     void jacobiIteration(const SP<ccsf> &b, const SP<const ccsf> &x) const
     {
	 jacobiIteration(*b, *x);
     }

     template<class FT>
     void jacobiIteration(FT &b, const FT &x) const;

     template<class FT>
     void jacobiIteration(const SP<FT> &b, const SP<const FT> &x) const
     {
	 jacobiIteration(*b, *x);
     }

     const ccsf &diagonal() const
     {
	 return diagonal_m;
     }

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_PCGDiffusionSolver

#endif                          // __PCGDiffusionSolver_MatrixP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of PCGDiffusionSolver/MatrixP1Diff.hh
//---------------------------------------------------------------------------//
