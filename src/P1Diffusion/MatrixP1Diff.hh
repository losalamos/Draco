//----------------------------------*-C++-*----------------------------------//
// MatrixP1Diff.hh
// Randy M. Roberts
// Mon Sep 28 08:33:48 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_MatrixP1Diff_hh__
#define __P1Diffusion_MatrixP1Diff_hh__

#include "ds++/SP.hh"

namespace rtt_P1Diffusion
{

 // Since rtt_P1Diffusion is a limited namespace
 // I risk putting in this using declaration for dsxx::SP.
 
 using dsxx::SP;
 
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
     
     SP<MT> spMesh;
     SP<const ccsf> spDiagonal;
     SP<const fcdsf> spOffDiagonal;
     FieldConstructor FC;

     // Mutable temporary data

     mutable SP<fcdsf> spxFC;
     mutable SP<fcdsf> spxSwap;
     mutable SP<fcdsf> spbFC;
     mutable SP<ccsf> spbtmp;
     mutable SP<ccsf> spxtmp;
     
   public:

     // CREATORS
    
     MatrixP1Diff(const SP<MT> &spMesh_,
		  const SP<const ccsf> &spDiagonal_,
		  const SP<const fcdsf> &spOffDiagonal_,
                  const FieldConstructor &FC_)
	 : spMesh(spMesh_), spDiagonal(spDiagonal_),
	   spOffDiagonal(spOffDiagonal_), FC(FC_)
     {
	 spxFC = new fcdsf(FC);
	 spxSwap = new fcdsf(FC);
	 spbFC = new fcdsf(FC);
	 spbtmp = new ccsf(FC);
	 spxtmp = new ccsf(FC);
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
	 multiply(*b, *x);
     }

     template<class FT>
     void jacobiIteration(FT &b, const FT &x) const;

     template<class FT>
     void jacobiIteration(const SP<FT> &b, const SP<const FT> &x) const
     {
	 multiply(*b, *x);
     }

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_MatrixP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/MatrixP1Diff.hh
//---------------------------------------------------------------------------//
