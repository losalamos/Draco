//----------------------------------*-C++-*----------------------------------//
// MatVecP1Diff.hh
// Randy M. Roberts
// Fri Sep 25 13:43:12 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __PCGDiffusionSolver_MatVecP1Diff_hh__
#define __PCGDiffusionSolver_MatVecP1Diff_hh__

#include "linalg/PCG_MatVec.hh"
#include "ds++/Mat.hh"
#include "ds++/SP.hh"

namespace rtt_PCGDiffusionSolver
{
 
 // Since rtt_PCGDiffusionSolver is a limited namespace
 // I risk putting in this using declaration for dsxx::SP, and dsxx::Mat1.
 
 using dsxx::SP;
 using dsxx::Mat1;

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
 class MatVecP1Diff : public PCG_MatVec<typename Matrix::value_type>
 {

     // NESTED CLASSES AND TYPEDEFS

   public:
     
     typedef typename Matrix::value_type value_type;
     typedef Mat1<value_type> Mat1;
     
     // DATA

   private:

     SP<const Matrix> spMatrix;
     
   public:

     // CREATORS
    
     MatVecP1Diff(const SP<const Matrix> &spMatrix_)
	 : spMatrix(spMatrix_)
     {
	 // empty
     }
     
     // MANIPULATORS
    
     void MatVec(Mat1 &b, const Mat1 &x);
     
     // ACCESSORS

     // none

   private:
    
     // IMPLEMENTATION

     // none
 };

} // namespace rtt_PCGDiffusionSolver

#endif                          // __PCGDiffusionSolver_MatVecP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of PCGDiffusionSolver/MatVecP1Diff.hh
//---------------------------------------------------------------------------//
