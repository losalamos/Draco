//----------------------------------*-C++-*----------------------------------//
// testSolverP1Diff.t.cc
// Randy M. Roberts
// Wed Sep 30 10:58:03 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "testSolverP1Diff.hh"
#include "../SolverP1Diff.hh"
#include "ds++/SP.hh"

namespace rtt_P1Diffusion_test
{

 using dsxx::SP;

 std::ostream &operator<<(std::ostream &os, const ccsf &rhs)
 {
     using std::endl;
    
     typedef ccsf FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
    
#if 0
     int iline = 0;
     for (FT::const_iterator it = rhs.begin(); it != rhs.end(); it++)
     {
	 os << std::setw(16) << *it << " ";
	 if (++iline % 6 == 0)
	     os << endl;
     }
     if (iline % 6 != 0)
	 os << endl;
#else
     os << endl;
     int icell = 0;
     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nz = rhs.get_Mesh().get_ncz();
     for (int k=0; k<nz; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 os << std::setw(4) << icell++ << ":";
		 os << " " << std::setw(16) << rhs(i,j,k);
		 os << endl;
	     }
#endif    

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 template<class MT>
 testSolverP1Diff<MT>::testSolverP1Diff(const Mesh_DB &mdb_,
					const pcg_DB &pcg_db_)
     : mdb(mdb_), pcg_db(pcg_db_)
 {
     // empty
 }
 
 template<class MT>
 void testSolverP1Diff<MT>::run()
 {
     SP<MT> spMesh(new MT(mdb));
     Solver solver(spMesh, pcg_db);

     SP<ccsf> spDiag(new ccsf(spMesh));
     *spDiag = -6.0;
    
     SP<fcdsf> spOffDiag(new fcdsf(spMesh));
     *spOffDiag = 1.0;

     SP<Matrix> spMatrix(new Matrix(spMesh, spDiag, spOffDiag));
    
     ccsf phi(spMesh);
     ccsf brhs(spMesh);
     ccsf phi0(spMesh);

     for (int i=0; i<brhs.size(); i++)
	 phi0[i] = i;

     spMatrix->multiply(brhs, phi0);
    
     solver.solve(phi, spMatrix, brhs);

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
     
     // std::cout << std::endl << phi << std::endl;
 }    

} // end namespace

//---------------------------------------------------------------------------//
//                              end of testSolverP1Diff.t.cc
//---------------------------------------------------------------------------//
