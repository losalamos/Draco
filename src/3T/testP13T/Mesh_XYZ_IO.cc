//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ_IO.cc
// Randy M. Roberts
// Fri Oct 16 14:15:14 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/Mesh_XYZ_IO.hh"
#include <iostream>
#include <iomanip>

namespace rtt_3T_testP13T
{

 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::ccsf &rhs)
 {
     using std::endl;
     
     typedef Mesh_XYZ::ccsf FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
     os << endl;
     
     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nz = rhs.get_Mesh().get_ncz();
     
     int icell = 0;
     for (int k=0; k<nz; k++)
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
     using std::endl;

     typedef Mesh_XYZ::cctf<std::vector<double> > FT;

     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);
    
     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nz = rhs.get_Mesh().get_ncz();
     
     int icell = 0;
     for (int k=0; k<nz; k++)
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
     using std::endl;

     typedef Mesh_XYZ::fcdsf FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nz = rhs.get_Mesh().get_ncz();
     
     int icell = 0;
     for (int k=0; k<nz; k++)
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
     using std::endl;

     typedef Mesh_XYZ::fcdtf<std::vector<double> > FT;
    
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::scientific << std::setprecision(6);

     os << endl;

     const int nx = rhs.get_Mesh().get_ncx();
     const int ny = rhs.get_Mesh().get_ncy();
     const int nz = rhs.get_Mesh().get_ncz();
     
     int icell = 0;
     for (int k=0; k<nz; k++)
	 for (int j=0; j<ny; j++)
	     for (int i=0; i<nx; i++)
	     {
		 for (int f=0; f<6; f++)
		 {
		     os << std::setw(5) << icell++ << ":" << f << ":";
		     int nmat = rhs(i,j,k,f).size();
		     for (int imat=0; imat<nmat; imat++)
			 os << " " << std::setw(16) << rhs(i,j,k,f)[imat];
		     os << endl;
		 }
	     }

     // restore the original flags.
    
     os.flags(fmtflags);

     return os;
 }

} // end namespace

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ_IO.cc
//---------------------------------------------------------------------------//
