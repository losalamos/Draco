//----------------------------------*-C++-*----------------------------------//
// GmvDump.t.cc
// Randy M. Roberts
// Wed Oct 14 10:19:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/GmvDump.hh"

#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Mat.hh"

#include <iostream>
#include <fstream>
#include <iomanip>

namespace rtt_3T_testP13T
{

 template<class MT>
 GmvDump<MT>::GmvDump(const std::string &fname_, const SP<MT> &spMesh_,
		       int cycle_, double time_)
     : fname(fname_), spMesh(spMesh_), nx(spMesh->get_ncx()),
       ny(spMesh->get_ncy()), nz(spMesh->get_ncz()), nzp(spMesh->get_nczp()),
       zoff(spMesh->get_zoff()),
       variablePrinted(false), cycle(cycle_), time(time_)
 {
     if (C4::node() != 0)
	 return;

     std::ofstream os(fname.c_str());
     
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

     VertId vid(nx, ny, nz);

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
     if (C4::node() != 0)
	 return;

     using std::endl;
    
     std::ofstream os(fname.c_str(), std::ios_base::app);

     if (variablePrinted)
	 os << "endvars" << endl;
     os << "probtime " << time << endl;
     os << "cycleno " << cycle << endl;
     os << "endgmv" << endl;
 }

 template<class MT>
 void GmvDump<MT>::dump(const typename MT::ccsf &var, const std::string &name)
{
     C4::HTSyncSpinLock htlock;

     using std::endl;

     std::ofstream os(fname.c_str(), std::ios_base::app);
     
     std::ios_base::fmtflags fmtflags = os.flags();

     os << std::setw(16) << std::scientific << std::setprecision(6);
	 
     if (C4::node() == 0 && !variablePrinted)
     {
	 os << "variable" << endl;
	 variablePrinted = true;
     }

     if (C4::node() == 0)
	 os << name << " 0" << endl;

     const int nitemsPerLine = 5;
     int nitems = 0;
     for (int k=zoff; k<zoff+nzp; k++)
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
