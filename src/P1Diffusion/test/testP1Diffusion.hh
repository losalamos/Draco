//----------------------------------*-C++-*----------------------------------//
// testP1Diffusion.hh
// Randy M. Roberts
// Thu Sep 24 16:45:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_test_testP1Diffusion_hh__
#define __P1Diffusion_test_testP1Diffusion_hh__

#include "P1Diffusion/P1Diffusion.hh"
#include "P1Diffusion/MatrixP1Diff.hh"
#include "P1Diffusion/MatVecP1Diff.hh"
#include "P1Diffusion/PreCondP1Diff.hh"
#include "P1Diffusion/SolverP1Diff.hh"

namespace rtt_P1Diffusion_test
{
 
 template<class MT>
 class testP1Diffusion
 {
     typedef rtt_P1Diffusion::SolverP1Diff<MT> MS;
     typedef MS MatrixSolver;
     typedef rtt_P1Diffusion::P1Diffusion<MT, MS> DiffSolver;

     typedef typename MT::ccsf ccsf;
     typedef typename MT::fcdsf fcdsf;
     typedef typename MT::bssf bssf;

     double D;
     double sigma;
     double q;
     double alphaBot;
     double alphaTop;
     double betaBot;
     double betaTop;
     double fBot;
     double fTop;
     double zBot;
     double zTop;
     double gamma;
     double A;
     double B;

     dsxx::SP<MT> spMesh;
     Diffusion_DB diffdb;
     pcg_DB pcg_db;

     int nx;
     int ny;
     int nz;
     
   public:

     testP1Diffusion(const dsxx::SP<MT> &spMesh_,
		     double D_, double sigma_, double q_,
		     double fTop_, double fBot_,
		     const Diffusion_DB &diffdb_,
		     const pcg_DB &pcg_db_);
     
     void run();

   private:

     void getAnalyticParams();
     void getSource(ccsf &Qsrc) const;
     void getPhi(ccsf &phi) const;
     double getPhi(double z) const;
     void setBoundary(bssf &alpha, bssf &beta, bssf &fb) const;
 };

} // end namespace

#endif                          // __P1Diffusion_test_testP1Diffusion_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/test/testP1Diffusion.hh
//---------------------------------------------------------------------------//
