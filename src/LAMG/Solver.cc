//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/Solver.cc
 * \author Randy M. Roberts
 * \date   Tue Jan 25 14:05:44 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Solver.hh"
#include "CompressedRowStorage.hh"
#include "LAMG_F90.hh"

namespace rtt_LAMG
{

void Solver::solve(DMat1 &x, const CompressedRowStorage &crs, const DMat1 &b)
{
    const IMat1 &rowPointer = crs.rowPointer();
    const IMat1 &colIndex = crs.colIndex();
    const DMat1 &val = crs.val();

    int nrows = rowPointer.size() - 1;
    int nentries = colIndex.size();

    Require(x.size() == nrows);
    Require(b.size() == nrows);

    // Error code, set to -1 for grins.
    // ierr == 0 --> no error
    // ierr  < 0 --> error
    // ierr  > 0 --> warning

    RTT_F90_INTEGER ierr = -1;

    LAMG_SOLVER_SOLVE(opaquePointer(), x.begin(), nrows, nentries,
		      b.begin(), rowPointer.begin(), colIndex.begin(),
		      val.begin(), ierr);
    Check(ierr >=0);
}

Solver::Representation::Representation(const Options &opts)
{
    IMat1 iopts(2);
    iopts[0] = opts.itsmax();
    iopts[1] = opts.levout();
    
    DMat1 dopts(1);
    dopts[0] = opts.tol();
    
    // Error code, set to -1 for grins.
    // ierr == 0 --> no error
    // ierr  < 0 --> error
    // ierr  > 0 --> warning

    RTT_F90_INTEGER ierr = -1;
    LAMG_SOLVER_CONSTRUCT(ptr, iopts.begin(), iopts.size(),
			  dopts.begin(), dopts.size(), ierr);
    Check(ierr >= 0);
}

Solver::Representation::~Representation()
{
    // Error code, set to -1 for grins.
    // ierr == 0 --> no error
    // ierr  < 0 --> error
    // ierr  > 0 --> warning

    RTT_F90_INTEGER ierr = -1;
    LAMG_SOLVER_DESTRUCT(ptr, ierr);
    Check(ierr >= 0);
}

} // end namespace rtt_LAMG

//---------------------------------------------------------------------------//
//                              end of Solver.cc
//---------------------------------------------------------------------------//
