//----------------------------------*-C++-*----------------------------------//
// PCG_Subroutines.cc
// Dave Nystrom
// Fri May  2 11:02:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "PCG_Subroutines.hh"

using namespace rtt_pcgWrap;

//---------------------------------------------------------------------------//
// xdfalt
//---------------------------------------------------------------------------//

void
rtt_pcgWrap::
xdfalt( int *iparm,
	float *fparm )
{
    sdfalt( iparm, fparm );
}

void
rtt_pcgWrap::
xdfalt( int *iparm,
	double *fparm )
{
    ddfalt( iparm, fparm );
}

//---------------------------------------------------------------------------//
// xbasr - basic iterative method
//---------------------------------------------------------------------------//

void
rtt_pcgWrap::
xbasr ( int& ijob,
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
	int& ier )
{
    sbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}

void
rtt_pcgWrap::
xbasr ( int& ijob,
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
	int& ier )
{
    dbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}

//---------------------------------------------------------------------------//
// xgmrsr - restarted gmres
//---------------------------------------------------------------------------//

void
rtt_pcgWrap::
xgmrsr( int& ijob,
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
	int& ier )
{
    sgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}

void
rtt_pcgWrap::
xgmrsr( int& ijob,
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
	int& ier )
{
    dgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}

//---------------------------------------------------------------------------//
// xcgr - conjugate gradient
//---------------------------------------------------------------------------//

void
rtt_pcgWrap::
xcgr( int& ijob,
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
      int& ier )
{
    scgr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	  ier );
}

void
rtt_pcgWrap::
xcgr( int& ijob,
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
      int& ier )
{
    dcgr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	  ier );
}

//---------------------------------------------------------------------------//
//                              end of PCG_Subroutines.cc
//---------------------------------------------------------------------------//
