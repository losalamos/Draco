//----------------------------------*-C++-*----------------------------------//
// PreCondP1Diff.hh
// Randy M. Roberts
// Mon Sep 28 16:30:41 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __PCGDiffusionSolver_PreCondP1Diff_hh__
#define __PCGDiffusionSolver_PreCondP1Diff_hh__

#include "linalg/PCG_PreCond.hh"
#include "ds++/Mat.hh"
#include "ds++/SP.hh"

namespace rtt_PCGDiffusionSolver
{

 // Since rtt_PCGDiffusionSolver is a limited namespace
 // I risk putting in this using declaration for dsxx::SP, and dsxx::Mat1.
 
 using dsxx::SP;
 using dsxx::Mat1;

 //===========================================================================//
 // class PreCondP1Diff - 
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
 class PreCondP1Diff : public PCG_PreCond<typename Matrix::value_type>
 {

     // NESTED CLASSES AND TYPEDEFS

   public:
     
     typedef typename Matrix::value_type value_type;
     typedef Mat1<value_type> Mat1;

     // DATA

     SP<const Matrix> spMatrix;
     int method;
    
   public:

     // CREATORS
    
     PreCondP1Diff(const SP<const Matrix> &spMatrix_, int method_)
	 : spMatrix(spMatrix_), method(method_)
     {
	 // empty
     }

     // MANIPULATORS
    
     void  Left_PreCond(Mat1 &x, Mat1 &b);
     void Right_PreCond(Mat1 &x, Mat1 &b);

     // ACCESSORS

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_PCGDiffusionSolver

#endif                          // __PCGDiffusionSolver_PreCondP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of PCGDiffusionSolver/PreCondP1Diff.hh
//---------------------------------------------------------------------------//
