//----------------------------------*-C++-*----------------------------------//
// testP1Diffusion.t.cc
// Randy M. Roberts
// Thu Sep 24 16:45:25 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "P1Diffusion/test/testP1Diffusion.hh"
#include <fstream>
#include <cmath>

namespace
{
 void solve2x2(double x[2], const double a[2][2], const double b[2])
 {
     const double det = a[0][0]*a[1][1] - a[0][1]*a[1][0];

     double ainv[2][2];

     ainv[0][0] = a[1][1] / det;
     ainv[1][1] = a[0][0] / det;
     ainv[0][1] = -a[0][1] / det;
     ainv[1][0] = -a[1][0] / det;

     for (int i=0; i<1; i++)
     {
	 x[i] = 0;
	 for (int j=0; j<1; j++)
	     x[i] += ainv[i][j]*b[j];
     }
 }
}

namespace rtt_P1Diffusion_test
{

 template<class MT>
 testP1Diffusion<MT>::testP1Diffusion(const dsxx::SP<MT> &spMesh_,
				      double D_, double sigma_, double q_,
				      double fTop_, double fBot_,
				      const Diffusion_DB &diffdb_,
				      const pcg_DB &pcg_db_)
     : spMesh(spMesh_), D(D_), sigma(sigma_), q(q_), fTop(fTop_), fBot(fBot_),
       diffdb(diffdb_), pcg_db(pcg_db_)
 {
     nx = spMesh->get_ncx();
     ny = spMesh->get_ncy();
     nz = spMesh->get_ncz();

     fcdsf zloc(spMesh);
     spMesh->get_zloc(zloc);
	 
     zBot = zloc(0,0,0,4);
     zTop = zloc(0,0,nz-1,5);

     alphaBot = diffdb.alpha_bottom;
     alphaTop = diffdb.alpha_top;
     betaBot = diffdb.beta_bottom;
     betaTop = diffdb.beta_top;

     std::cout << " alphaBot, betaBot: " 
	       << alphaBot << " " << betaBot
	       << " alphaTop, betaTop: "
	       << alphaTop << " " << betaTop << std::endl;
     
     getAnalyticParams();
 }

 template<class MT>
 void testP1Diffusion<MT>::run()
 {
     SP<MatrixSolver> spMatrixSolver = new MatrixSolver(spMesh, pcg_db);
     SP<DiffSolver> spDiffSolver = new DiffSolver(diffdb, spMesh,
						  spMatrixSolver);

     fcdsf DFC(spMesh);
     DFC = D;

     ccsf sigmaCC(spMesh);
     sigmaCC = sigma;

     ccsf QCC(spMesh);
     getSource(QCC);

     fcdsf Fprime(spMesh);
     Fprime = 0.0;

     bssf alpha(spMesh);
     bssf beta(spMesh);
     bssf f(spMesh);

     setBoundary(alpha, beta, f);
     
     ccsf phi(spMesh);
     fcdsf F(spMesh);

     spDiffSolver->solve(phi, F, DFC, sigmaCC, QCC, Fprime, alpha, beta, f);

     ccsf phi0(spMesh);
     getPhi(phi0);

     double error = 0.0;
     double l2nPhi0 = 0.0;
     for (int i=0; i<phi.size(); i++)
     {
	 error += (phi[i]-phi0[i])*(phi[i]-phi0[i]);
	 l2nPhi0 += phi0[i]*phi0[i];
     }

     error /= l2nPhi0;
     error = std::sqrt(error);

     std::cout << "error: " << error << std::endl;

     std::cout << "phi(0), F(0): " << phi[0] << " " << F(0,0,0,4)
	       << std::endl;
     std::cout << "phi(nz-1), F(nz-1): " << phi[nz-1] << " " << F(0,0,nz-1,5)
	       << std::endl;

     ccsf zc(spMesh);
     spMesh->get_zloc(zc);
     std::ofstream ofs("testP1Diffusion.dat");
     for (int i=0; i<phi.size(); i++)
	 ofs << zc[i] << "\t" << phi0[i] << "\t" << phi[i] << std::endl;
 }
 
 template<class MT>
 void testP1Diffusion<MT>::getAnalyticParams()
 {
     gamma = std::sqrt(sigma/D);

#if 0
     double G = 4.*J0 - 2.*D*Q/(sigma*sigma);
     double H = 4.*J1 - Q/sigma*(2.*D*(2.+1./sigma)+1.);
  
     double A11 = (1.-2.*D*gamma);
     double A12 = (1.+2.*D*gamma);
     double A21 = exp(gamma)*A12;
     double A22 = exp(-gamma)*A11;

     double Det = (A11*A22 - A12*A21);

     double A = (A22*G - A12*H)/ Det;
     double B = (-A21*G + A11*H)/ Det;
#endif
     
     const double rhsBot = fBot - alphaBot*q/sigma*(zBot*zBot + 2.0*D/sigma)
	 - betaBot*D*2.0*q*zBot/sigma;

     const double expGammaBot = std::exp(gamma*zBot);
     const double ACoefBot = ( betaBot*D*gamma + alphaBot) * expGammaBot;
     const double BCoefBot = (-betaBot*D*gamma + alphaBot) / expGammaBot;

     const double rhsTop = fTop - alphaTop*q/sigma*(zTop*zTop + 2.0*D/sigma)
	 + betaTop*D*2.0*q*zTop/sigma;

     const double expGammaTop = std::exp(gamma*zTop);
     const double ACoefTop = (-betaTop*D*gamma + alphaTop) * expGammaTop;
     const double BCoefTop = ( betaTop*D*gamma + alphaTop) / expGammaTop;

     double A11 = ACoefBot;
     double A12 = BCoefBot;
     double A21 = ACoefTop;
     double A22 = BCoefTop;

     double G = rhsBot;
     double H = rhsTop;

     const double Det = (A11*A22 - A12*A21);

     A = (A22*G - A12*H)/ Det;
     B = (-A21*G + A11*H)/ Det;

     std::cout << "A, B: " << A << " " << B << std::endl;
 }
 
 template<class MT>
 void testP1Diffusion<MT>::getSource(ccsf &Qsrc) const
 {
     ccsf zc(spMesh);
     spMesh->get_zloc(zc);

     for (int i=0; i<Qsrc.size(); i++)
     {
	 const double r2 = zc[i]*zc[i];
	 Qsrc[i] = q*r2;
     }
 }

 template<class MT>
 void testP1Diffusion<MT>::getPhi(ccsf &phi) const
 {
     ccsf zc(spMesh);
     spMesh->get_zloc(zc);

     for (int i=0; i<phi.size(); i++)
     {
	 phi[i] = getPhi(zc[i]);
     }
 }

 template<class MT>
 double testP1Diffusion<MT>::getPhi(double z) const
 {
     const double expGamma = std::exp(gamma*z);

     double phi = A*expGamma + B/expGamma + q*z*z/sigma + 2.0*D*q/(sigma*sigma);
     return phi;
 }

 template<class MT>
 void testP1Diffusion<MT>::setBoundary(bssf &alpha, bssf &beta,
				       bssf &fb) const
 {
     alpha = 0.0;
     beta = 1.0;
     fb = 0.0;
     
     for (int i=0; i<nx; i++)
     {
	 for (int j=0; j<ny; j++)
	 {
	     alpha(i, j, 0   , 4) = alphaBot;
	     beta (i, j, 0   , 4) = betaBot;
	     fb   (i, j, 0   , 4) = fBot;
	     alpha(i, j, nz-1, 5) = alphaTop;
	     beta (i, j, nz-1, 5) = betaTop;
	     fb   (i, j, nz-1, 5) = fTop;
	 }
     }
 }    

} // end namespace

//---------------------------------------------------------------------------//
//                              end of testP1Diffusion.t.cc
//---------------------------------------------------------------------------//
