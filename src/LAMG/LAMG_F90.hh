//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/LAMG_F90.hh
 * \author Randy M. Roberts
 * \date   Tue Jan 25 14:15:05 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_LAMG_F90_hh__
#define __LAMG_LAMG_F90_hh__

#include <LAMG/config.h>
#include "c2f.h"

#define LAMG_SOLVER_CONSTRUCT(ptr, iopts, niopts, dopts, ndopts, ierr) \
      CCALLSFSUB6(LAMG_SOLVER_CONSTRUCT, lamg_solver_construct, \
                  INT, ptr, INT_PTR, iopts, INT, niopts, \
                  DBL_PTR, dopts, INT, ndopts, INT, ierr)

#define LAMG_SOLVER_DESTRUCT(ptr, ierr) \
      CCALLSFSUB2(LAMG_SOLVER_DESTRUCT, lamg_solver_destruct, \
                  INT, ptr, INT, ierr)

#define LAMG_SOLVER_SOLVE(ptr, x, nrows, nentries, b, rowPointer, \
                          colIndex, val, ierr) \
      CCALLSFSUB9(LAMG_SOLVER_SOLVE, lamg_solver_solve, \
                  INT, ptr, DBL_PTR, x, INT, nrows, INT, nentries, \
                  DBL_PTR, b, INT_PTR, rowPointer, INT_PTR, colIndex, \
                  DBL_PTR, val, INT, ierr)

namespace rtt_LAMG
{

extern "C" {
    void LAMG_SOLVER_CONSTRUCT(RTT_F90_INTEGER &ptr,
			       const RTT_F90_INTEGER *iopts,
			       const RTT_F90_INTEGER &niopts,
			       const RTT_F90_DOUBLE *dopts,
			       const RTT_F90_INTEGER &ndopts,
			       RTT_F90_INTEGER &ierr);
    
    void LAMG_SOLVER_DESTRUCT(const RTT_F90_INTEGER &ptr,
			      RTT_F90_INTEGER &ierr);
    
    void LAMG_SOLVER_SOLVE(const RTT_F90_INTEGER &ptr,
			   RTT_F90_DOUBLE *x,
			   const RTT_F90_INTEGER &nrows,
			   const RTT_F90_INTEGER &nentries,
			   const RTT_F90_DOUBLE *b,
			   const RTT_F90_INTEGER *rowPointer,
			   const RTT_F90_INTEGER *colIndex,
			   const RTT_F90_DOUBLE *val,
			   RTT_F90_INTEGER &ierr);
}

} // end namespace rtt_LAMG

#endif                          // __LAMG_LAMG_F90_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/LAMG_F90.hh
//---------------------------------------------------------------------------//
