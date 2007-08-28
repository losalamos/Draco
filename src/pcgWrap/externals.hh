//----------------------------------*-C++-*----------------------------------//
// externals.hh
// Rob Lowrie
// Tue Aug 28 06:45:57 2007
//---------------------------------------------------------------------------//

// Defines actual prototypes for routines in the external PCG library.

#ifndef __pcgWrap_externals_hh__
#define __pcgWrap_externals_hh__

#include <pcgWrap/config.h>

// Handle machine specific Fortran name linkage.

#if defined(sun) || defined(__sun) || defined(__sgi) || defined(__linux) || defined(__alpha) || defined(__GNUC__)

#define sdfalt sdfalt_
#define ddfalt ddfalt_
#define cdfalt cdfalt_
#define zdfalt zdfalt_

#define sbasr  sbasr_
#define dbasr  dbasr_
#define cbasr  cbasr_
#define zbasr  zbasr_

#define sgmrsr sgmrsr_
#define dgmrsr dgmrsr_
#define cgmrsr cgmrsr_
#define zgmrsr zgmrsr_

#define scgr   scgr_
#define dcgr   dcgr_
#define ccgr   ccgr_
#define zcgr   zcgr_
#endif

// Function prototypes for PCG f77 subroutines.

#ifdef USE_PCGLIB

extern "C" {

// ... then PCGLIB is actually supported

    void sdfalt( int *iparm,
		 float *fparm );
    
    void ddfalt( int *iparm,
		 double *fparm );

// Basic iterative method.
    void sbasr ( int& ijob,
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
    
    void dbasr ( int& ijob,
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
    void sgmrsr( int& ijob,
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
    
    void dgmrsr( int& ijob,
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
    void scgr  ( int& ijob,
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
    
    void dcgr  ( int& ijob,
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

} // extern C

#else  // ... then PCGLIB is not supported

    // Make these functions no-ops
    void sdfalt( int *iparm,
		 float *fparm ) {}
    
    void ddfalt( int *iparm,
		 double *fparm ) {}

// Basic iterative method.
    void sbasr ( int& ijob,
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
		 int& ier ) {}
    
    void dbasr ( int& ijob,
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
		 int& ier ) {}

// Restarted GMRES.
    void sgmrsr( int& ijob,
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
		 int& ier ) {}
    
    void dgmrsr( int& ijob,
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
		 int& ier ) {}

// Conjugate gradient.
    void scgr  ( int& ijob,
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
		 int& ier ) {}
    
    void dcgr  ( int& ijob,
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
		 int& ier ) {}
#endif // USE_PCGLIB

#endif                          // __pcgWrap_externals_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/externals.hh
//---------------------------------------------------------------------------//
