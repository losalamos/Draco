//----------------------------------*-C++-*----------------------------------//
// PCG_Subroutines.hh
// Dave Nystrom
// Fri May  2 11:02:51 1997
//---------------------------------------------------------------------------//
// @> Defines wrapped PCG routine prototypes.
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_PCG_Subroutines_hh__
#define __pcgWrap_PCG_Subroutines_hh__

namespace rtt_pcgWrap {

void xdfalt( int *iparm,
	     float *fparm );

void xdfalt( int *iparm,
	     double *fparm );

// Basic iterative method.
void xbasr ( int& ijob,
	     int& ireq,
	     float *x,
	     float *xex,
	     const float *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     float *fwk,
	     int *iparm,
	     float *fparm,
	     int& ier );

void xbasr ( int& ijob,
	     int& ireq,
	     double *x,
	     double *xex,
	     const double *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     double *fwk,
	     int *iparm,
	     double *fparm,
	     int& ier );

// Restarted GMRES.
void xgmrsr( int& ijob,
	     int& ireq,
	     float *x,
	     float *xex,
	     const float *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     float *fwk,
	     int *iparm,
	     float *fparm,
	     int& ier );

void xgmrsr( int& ijob,
	     int& ireq,
	     double *x,
	     double *xex,
	     const double *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     double *fwk,
	     int *iparm,
	     double *fparm,
	     int& ier );

// Conjugate gradient.
void xcgr  ( int& ijob,
	     int& ireq,
	     float *x,
	     float *xex,
	     const float *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     float *fwk,
	     int *iparm,
	     float *fparm,
	     int& ier );

void xcgr  ( int& ijob,
	     int& ireq,
	     double *x,
	     double *xex,
	     const double *b,
	     int& iva,
	     int& ivql,
	     int& ivqr,
	     int *iwk,
	     double *fwk,
	     int *iparm,
	     double *fparm,
	     int& ier );

} // namespace rtt_pcgWrap

#endif                          // __pcgWrap_PCG_Subroutines_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/PCG_Subroutines.hh
//---------------------------------------------------------------------------//
