//----------------------------------*-C++-*----------------------------------//
// SolverP1Diff.hh
// Randy M. Roberts
// Tue Sep 29 16:11:17 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __PCGDiffusionSolver_SolverP1Diff_hh__
#define __PCGDiffusionSolver_SolverP1Diff_hh__

#include "pcg_DB.hh"
#include "ds++/SP.hh"
#include "linalg/PCG_Ctrl.hh"
#include "MatrixP1Diff.hh"

namespace rtt_PCGDiffusionSolver
{

 // Within namespace forward declarations.
 
 // template<class MT> class MatrixP1Diff;
 template<class Matrix> class MatVecP1Diff;
 template<class Matrix> class PreCondP1Diff;

 // Inside a small namespace I feel I can put in a using declaration.
 
 using rtt_dsxx::SP;
 
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
     typedef typename MT::FieldConstructor FieldConstructor;
     typedef MatrixP1Diff<MT> Matrix;
     typedef MatVecP1Diff<Matrix> MatVec;
     typedef PreCondP1Diff<Matrix> PreCond;
     typedef typename ccsf::value_type value_type;

   private:
     
     typedef PCG_Ctrl<value_type> PCG_Ctrl;

     // DATA

   private:
     
     SP<const MT> spMesh;
     PCG_Ctrl pcg_ctrl;
     // precond method 0 - none, 1 - left, 2 - right, 3 - both
     int iqside;
    
   public:

     // CREATORS
    
     SolverP1Diff(const SP<const MT>& spMesh, const FieldConstructor &fCtor_in,
                  const pcg_DB& pcg_db);

     // MANIPULATORS

     void solve(ccsf &phi, const SP<const Matrix> &spMatrix,
		const ccsf &brhs);
     
     // ACCESSORS

   private:
    
     // IMPLEMENTATION

     // STATIC IMPLEMENTATION

     static typename PCG_Ctrl::Method pcgMethod(int itmeth);
     static typename PCG_Ctrl::OutputLevel outputLevel(int levout);
     static typename PCG_Ctrl::StopTest stopTest(int ntest);
     static typename PCG_Ctrl::Uinit uninit(int iuinit);
     static typename PCG_Ctrl::Precon precon(int iqside);
     static typename PCG_Ctrl::Logical logical(int value);
     
 };

} // end namespace rtt_PCGDiffusionSolver

#endif                          // __PCGDiffusionSolver_SolverP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of PCGDiffusionSolver/SolverP1Diff.hh
//---------------------------------------------------------------------------//
