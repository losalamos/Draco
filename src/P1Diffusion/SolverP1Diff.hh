//----------------------------------*-C++-*----------------------------------//
// SolverP1Diff.hh
// Randy M. Roberts
// Tue Sep 29 16:11:17 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_SolverP1Diff_hh__
#define __P1Diffusion_SolverP1Diff_hh__

#include "ds++/SP.hh"
#include "linalg/PCG_Ctrl.hh"
#include "P1Diffusion/MatrixP1Diff.hh"

namespace rtt_P1Diffusion
{

 // Within namespace forward declarations.
 
 // template<class MT> class MatrixP1Diff;
 template<class Matrix> class MatVecP1Diff;
 template<class Matrix> class PreCondP1Diff;

 // Inside a small namespace I feel I can put in a using declaration.
 
 using dsxx::SP;
 
 //===========================================================================//
 // class SolverP1Diff - 
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
 class SolverP1Diff
 {
     // NESTED CLASSES AND TYPEDEFS

   public:
     
     typedef typename MT::ccsf ccsf;
     typedef MatrixP1Diff<MT> Matrix;
     typedef MatVecP1Diff<Matrix> MatVec;
     typedef PreCondP1Diff<Matrix> PreCond;

     // DATA

   private:
     
     SP<const MT> spMesh;
     PCG_Ctrl<typename ccsf::value_type> pcg_ctrl;
    
   public:

     // CREATORS
    
     SolverP1Diff(const SP<const MT>& spMesh_, const pcg_DB& pcg_db);

     // MANIPULATORS

     void solve(ccsf &phi, const SP<const Matrix> &spMatrix,
		const ccsf &brhs);
     
     // ACCESSORS

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_SolverP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/SolverP1Diff.hh
//---------------------------------------------------------------------------//
