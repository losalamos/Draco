//----------------------------------*-C++-*----------------------------------//
// GmvDump.t.cc
// Randy M. Roberts
// Wed Oct 14 10:19:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/GmvDump.hh"

#include <iostream>
#include <iomanip>

namespace rtt_3T_testP13T
{

 template<class MT>
 GmvDump<MT>:: GmvDump(std::ostream &os_, const SP<MT> &spMesh_,
		       int cycle_, double time_)
     : os(os_), spMesh(spMesh_), nx(spMesh->get_ncx()),
       ny(spMesh->get_ncy()), nz(spMesh->get_ncz()), vid(nx, ny, nz),
       variablePrinted(false), cycle(cycle_), time(time_)
 {
     using std::endl;
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::setw(16) << std::scientific << std::setprecision(6);

     os << "gmvinput ascii" << endl;

     double dx = spMesh->get_dx();
     double dy = spMesh->get_dy();
     double dz = spMesh->get_dz();

     int nnodes = (nx+1)*(ny+1)*(nz+1);
     os << "nodes " << nnodes << endl;

     const int nitemsPerLine = 5;
     int nitems = 0;

     for (int k=0; k<nz+1; k++)
     {
	 for (int j=0; j<ny+1; j++)
	 {
	     for (int i=0; i<nx+1; i++)
	     {
		 os << i*dx;
		 if (++nitems % nitemsPerLine)
		     os << " ";
		 else
		     os << endl;
	     }
	 }
     }
     if (nitems % nitemsPerLine)
	 os << endl;

     nitems = 0;
     for (int k=0; k<nz+1; k++)
     {
	 for (int j=0; j<ny+1; j++)
	 {
	     for (int i=0; i<nx+1; i++)
	     {
		 os << j*dy;
		 if (++nitems % nitemsPerLine)
		     os << " ";
		 else
		     os << endl;
	     }
	 }
     }
     if (nitems % nitemsPerLine)
	 os << endl;

     nitems = 0;
     for (int k=0; k<nz+1; k++)
     {
	 for (int j=0; j<ny+1; j++)
	 {
	     for (int i=0; i<nx+1; i++)
	     {
		 os << k*dz;
		 if (++nitems % nitemsPerLine)
		     os << " ";
		 else
		     os << endl;
	     }
	 }
     }
     if (nitems % nitemsPerLine)
	 os << endl;

     os << "cells " << nx*ny*nz << endl;

     for (int k=0; k<nz; k++)
     {
	 for (int j=0; j<ny; j++)
	 {
	     for (int i=0; i<nx; i++)
	     {
		 os << "hex 8" << endl;
		 os << vid(i,   j,   k+1) << " ";
		 os << vid(i+1, j,   k+1) << " ";
		 os << vid(i+1, j+1, k+1) << " ";
		 os << vid(i,   j+1, k+1) << " ";
		 os << vid(i,   j,   k  ) << " ";
		 os << vid(i+1, j,   k  ) << " ";
		 os << vid(i+1, j+1, k  ) << " ";
		 os << vid(i,   j+1, k  ) << " ";
		 os << endl;
	     }
	 }
     }

     // restore the original flags.
	 
     os.flags(fmtflags);
 }

 template<class MT>
 GmvDump<MT>::~GmvDump()
 {
     using std::endl;
    
     if (variablePrinted)
	 os << "endvars" << endl;
     os << "probtime " << time << endl;
     os << "cycleno " << cycle << endl;
     os << "endgmv" << endl;
 }

 template<class MT>
 void GmvDump<MT>::dump(const typename MT::ccsf &var, const std::string &name)
 {
     using std::endl;
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::setw(16) << std::scientific << std::setprecision(6);
	 
     if (!variablePrinted)
     {
	 os << "variable" << endl;
	 variablePrinted = true;
     }
	 
     os << name << " 0" << endl;

     const int nitemsPerLine = 5;
     int nitems = 0;
     for (int k=0; k<nz; k++)
     {
	 for (int j=0; j<ny; j++)
	 {
	     for (int i=0; i<nx; i++)
	     {
		 os << var(i,j,k);
		 if (++nitems % nitemsPerLine)
		     os << " ";
		 else
		     os << endl;
	     }
	 }
     }

     if (nitems % nitemsPerLine)
	 os << endl;

	 // restore the original flags.
	 
     os.flags(fmtflags);
 }

} // end namespace rtt_3T_testP13T

//---------------------------------------------------------------------------//
//                              end of GmvDump.t.cc
//---------------------------------------------------------------------------//
