//----------------------------------*-C++-*----------------------------------//
// testSolverP1Diff.hh
// Randy M. Roberts
// Wed Sep 30 10:58:03 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_test_testSolverP1Diff_hh__
#define __P1Diffusion_test_testSolverP1Diff_hh__

#include "P1Diffusion/SolverP1Diff.hh"
#include "mesh/Mesh_DB.hh"
#include "linalg/pcg_DB.hh"

namespace rtt_P1Diffusion_test
{
 //===========================================================================//
 // class testSolverP1Diff - 
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
 class testSolverP1Diff
 {

     // NESTED CLASSES AND TYPEDEFS

     typedef typename MT::ccsf ccsf;
     typedef typename MT::fcdsf fcdsf;

     typedef rtt_P1Diffusion::SolverP1Diff<MT> Solver;
     typedef typename Solver::Matrix Matrix;

     // DATA

     Mesh_DB mdb;
     pcg_DB pcg_db;
     
   public:

     // CREATORS

     testSolverP1Diff(const Mesh_DB &mdb_, const pcg_DB &pcg_db_);
     
     // MANIPULATORS

     void run();
    
     // ACCESSORS

   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_P1Diffusion_test

#endif                          // __P1Diffusion_test_testSolverP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/test/testSolverP1Diff.hh
//---------------------------------------------------------------------------//
