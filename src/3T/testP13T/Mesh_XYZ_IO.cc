//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ_IO.cc
// Randy M. Roberts
// Fri Oct 16 14:15:14 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/Mesh_XYZ_IO.hh"
#include "3T/testP13T/utils.hh"
#include "ds++/Mat.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Mat.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

namespace rtt_3T_testP13T
{

 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::ccsf &rhs)
 {
     C4::HTSyncSpinLock h;
     
     using std::endl;
     
     typedef Mesh_XYZ::ccsf FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
     os << endl;
     
     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nzp = rhs.get_Mesh().get_nczp();
     const int zoff = rhs.get_Mesh().get_zoff();
     
     int icell = rhs.get_Mesh().get_goff();
     
     for (int k=zoff; k<zoff+nzp; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 os << std::setw(5) << icell++ << ":";
		 os << " " << std::setw(16) << rhs(i,j,k);
		 os << endl;
	     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os, 
			  const Mesh_XYZ::cctf<std::vector<double> > &rhs)
 {
     C4::HTSyncSpinLock h;

     using std::endl;

     typedef Mesh_XYZ::cctf<std::vector<double> > FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
    
     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nzp = rhs.get_Mesh().get_nczp();
     const int zoff = rhs.get_Mesh().get_zoff();
     
     int icell = rhs.get_Mesh().get_goff();
     for (int k=zoff; k<zoff+nzp; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 os << std::setw(5) << icell++ << ":";
		 int nmat = rhs(i,j,k).size();
		 for (int imat = 0; imat < nmat; imat++)
		 {
		     os << " " << std::setw(16) << rhs(i,j,k)[imat];
		 }
		 os << endl;
	     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::fcdsf &rhs)
 {
     C4::HTSyncSpinLock h;

     using std::endl;

     typedef Mesh_XYZ::fcdsf FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nzp = rhs.get_Mesh().get_nczp();
     const int zoff = rhs.get_Mesh().get_zoff();
     
     int icell = rhs.get_Mesh().get_goff();
     for (int k=zoff; k<zoff+nzp; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 os << std::setw(5) << icell++ << ":";
		 for (int f=0; f<6; f++)
		     os << " " << std::setw(16) << rhs(i,j,k,f);
		 os << endl;
	     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os,
			  const Mesh_XYZ::fcdtf<std::vector<double> > &rhs)
 {
     C4::HTSyncSpinLock h;

     using std::endl;

     typedef Mesh_XYZ::fcdtf<std::vector<double> > FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nzp = rhs.get_Mesh().get_nczp();
     const int zoff = rhs.get_Mesh().get_zoff();
     
     int icell = rhs.get_Mesh().get_goff();
     for (int k=zoff; k<zoff+nzp; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 for (int f=0; f<6; f++)
		 {
		     os << std::setw(5) << icell << ":" << f << ":";
		     int nmat = rhs(i,j,k,f).size();
		     for (int imat=0; imat<nmat; imat++)
			 os << " " << std::setw(16) << rhs(i,j,k,f)[imat];
		     os << endl;
		 }
		 icell++;
	     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::bssf &rhs)
 {
     C4::HTSyncSpinLock h;

     using std::endl;
     
     typedef Mesh_XYZ::bssf FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
     os << endl;
     
     int bndry = 0;
     for (FT::const_iterator bit = rhs.begin();
	  bit != rhs.end();
	  bit++)
     {
	 os << std::setw(5) << bndry++ << ":";
	 os << " " << std::setw(16) << *bit;
	 os << endl;
     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

 void setBoundary(Mesh_XYZ::bssf &bndry, double left, double right, double
		  front, double back, double bottom, double top)
 {
     const Mesh_XYZ &mesh = bndry.get_Mesh();
     
     const int nx = mesh.get_ncx();
     const int ny = mesh.get_ncy();
     const int nzp = mesh.get_nczp();
     const int zoff = mesh.get_zoff();
     
     for (int j=0; j<ny; j++)
     {
	 for (int k=zoff; k<zoff+nzp; k++)
	 {
	     bndry(0,    j, k, LEFT) = left;
	     bndry(nx-1, j, k, RIGHT) = right;
	 }
     }

     for (int i=0; i<nx; i++)
     {
	 for (int k=zoff; k<zoff+nzp; k++)
	 {
	     bndry(i, 0   , k, FRONT) = front;
	     bndry(i, ny-1, k, BACK) = back;
	 }
     }

     if (C4::node() == 0)
     {
	 for (int i=0; i<nx; i++)
	     for (int j=0; j<ny; j++)
		 bndry(i, j, 0   , BOTTOM) = bottom;
     }

     if (C4::node() == C4::nodes()-1)
     {
	 for (int i=0; i<nx; i++)
	     for (int j=0; j<ny; j++)
		 bndry(i, j, zoff+nzp-1, TOP) = top;
     }
 }

 void setBoundary(Mesh_XYZ::bstf<std::vector<double> > &bndry,
		  double val, Faces face)
 {
     // RMR this code is Mesh_XYZ specific!!!
    
     const Mesh_XYZ &mesh = bndry.get_Mesh();
     
     const int nx = mesh.get_ncx();
     const int ny = mesh.get_ncy();
     const int nzp = mesh.get_nczp();
     const int zoff = mesh.get_zoff();

     using std::fill;
     
     switch(face)
     {
     case LEFT:
	 for (int j=0; j<ny; j++)
	 {
	     for (int k=zoff; k<zoff+nzp; k++)
	     {
		 fill(bndry(0, j, k, LEFT).begin(),
		      bndry(0, j, k, LEFT).end(), val);
	     }
	 }
	 break;
     case RIGHT:
	 for (int j=0; j<ny; j++)
	 {
	     for (int k=zoff; k<zoff+nzp; k++)
	     {
		 fill(bndry(nx-1, j, k, RIGHT).begin(),
		      bndry(nx-1, j, k, RIGHT).end(), val);
	     }
	 }
	 break;
     case FRONT:
	 for (int i=0; i<nx; i++)
	 {
	     for (int k=zoff; k<zoff+nzp; k++)
	     {
		 fill(bndry(i, 0 , k, FRONT).begin(),
		      bndry(i, 0 , k, FRONT).end(), val);
	     }
	 }
	 break;
     case BACK:
	 for (int i=0; i<nx; i++)
	 {
	     for (int k=zoff; k<zoff+nzp; k++)
	     {
		 fill(bndry(i, ny-1, k, BACK).begin(),
		      bndry(i, ny-1, k, BACK).end(), val);
	     }
	 }
	 break;
     case BOTTOM:
	 if (C4::node() == 0)
	 {
	     for (int i=0; i<nx; i++)
		 for (int j=0; j<ny; j++)
		     fill(bndry(i, j, 0 , BOTTOM).begin(),
			  bndry(i, j, 0 , BOTTOM).end(), val);
	 }
	 break;
     case TOP:
	 if (C4::node() == C4::nodes()-1)
	 {
	     for (int i=0; i<nx; i++)
		 for (int j=0; j<ny; j++)
		     fill(bndry(i, j, zoff+nzp-1, TOP).begin(),
			  bndry(i, j, zoff+nzp-1, TOP).end(), val);
	 }
	 break;
     default:
	 Assert(false);
	 break;
     }
 }

 void setTempFromFile(Mesh_XYZ::ccsf &Temp, const std::string &filename,
		      double floor)
 {
     std::ifstream ifs(filename.c_str());

     ifs.ignore(10000, '\n');
     ifs.ignore(10000, '\n');

     const Mesh_XYZ &mesh = Temp.get_Mesh();
     
     const int nx = mesh.get_ncx();
     const int ny = mesh.get_ncy();
     const int nz = mesh.get_ncz();
     const int nzp = mesh.get_nczp();
     const int zoff = mesh.get_zoff();

     std::vector<double> T1(nz);
     
     for (int k=0; k<nz; k++)
     {
	 double r, T2, T3;

	 ifs >> r >> T1[k] >> T2 >> T3;
	 if (ifs.eof())
	     throw std::runtime_error("Premature EOF reading temperatures.");

	 if (T1[k] < floor)
	     T1[k] = floor;
     }

     for (int k=zoff; k<zoff+nzp; k++)
	 for (int i=0; i<nx; i++)
	     for (int j=0; j<ny; j++)
		 Temp(i, j, k) = T1[k];
 }

 void dumpInZ(const std::string &fname, int cycle, double time,
	      const Mesh_XYZ::ccsf &Temp)
 {
     C4::HTSyncSpinLock htlock;

     using std::endl;

     const Mesh_XYZ &mesh = Temp.get_Mesh();
     
     const int nzp = mesh.get_nczp();
     const int zoff = mesh.get_zoff();

     std::ios_base::openmode openmode = std::ios_base::app;
     if (C4::node() == 0)
	 openmode = std::ios_base::out | std::ios_base::trunc;
	    
     std::ofstream ofs(fname.c_str(), openmode);

     ofs << "# testFullP13T cycle=" << cycle << " time=" << time << endl;
     ofs << "# z \t TRad(z)" << endl;

     const double dz = mesh.get_dz();
     for (int m=zoff; m<nzp+zoff; m++)
	 ofs << dz*(m+.5) << '\t' << Temp(0,0,m)
	     << endl;
 }
 
} // end namespace

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ_IO.cc
//---------------------------------------------------------------------------//
