//----------------------------------*-C++-*----------------------------------//
// SolverP1Diff.cc
// Randy M. Roberts
// Tue Sep 29 16:11:17 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/SolverP1Diff.hh"
#include "P1Diffusion/MatrixP1Diff.hh"
#include "P1Diffusion/MatVecP1Diff.hh"
#include "P1Diffusion/PreCondP1Diff.hh"

namespace rtt_P1Diffusion
{

 template<class MT>
 SolverP1Diff<MT>::SolverP1Diff(const SP<const MT>& spMesh_,
				const pcg_DB& pcg_db)
     : spMesh(spMesh_), pcg_ctrl(pcg_db, spMesh_->get_ncells())
 {
     // empty
 }

 template<class MT>
 void SolverP1Diff<MT>::solve(ccsf &phi, const SP<const Matrix> &spMatrix,
			      const ccsf &brhs)
 {
     SP<MatVec> spMatVec = new MatVec(spMatrix);
     
     const int method = 1;
     SP<PreCond> spPreCond = new PreCond(spMatrix, method);

     // Now solve the matrix equation A.x = rhs.

     dsxx::Mat1<double> phitmp(phi.size());
     std::copy(phi.begin(), phi.end(), phitmp.begin());
    
     dsxx::Mat1<double> brhstmp(brhs.size());
     std::copy(brhs.begin(), brhs.end(), brhstmp.begin());

     pcg_ctrl.pcg_fe(phitmp, brhstmp, spMatVec, spPreCond);
     
     std::copy(phitmp.begin(), phitmp.end(), phi.begin());
 }

} // end namespace



//---------------------------------------------------------------------------//
//                              end of SolverP1Diff.cc
//---------------------------------------------------------------------------//
