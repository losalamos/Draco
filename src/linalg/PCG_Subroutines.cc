//----------------------------------*-C++-*----------------------------------//
// PCG_Subroutines.cc
// Dave Nystrom
// Fri May  2 11:02:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/PCG_Subroutines.hh"
using namespace pcg;

//---------------------------------------------------------------------------//
// xdfalt
//---------------------------------------------------------------------------//

void pcg::xdfalt( int *iparm, float *fparm ) {
    sdfalt( iparm, fparm );
}
void pcg::xdfalt( int *iparm, double *fparm ) {
    ddfalt( iparm, fparm );
}
void pcg::xdfalt( int *iparm, complex<float> *fparm ) {
    cdfalt( iparm, fparm );
}
void pcg::xdfalt( int *iparm, complex<double> *fparm ) {
    zdfalt( iparm, fparm );
}

//---------------------------------------------------------------------------//
// xbasr
//---------------------------------------------------------------------------//

void pcg::xbasr ( int& ijob, int& ireq,
		  float *x, float *xex, float *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  float *fwk, int *iparm, float *fparm,
		  int& ier ) {
    sbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}
void pcg::xbasr ( int& ijob, int& ireq,
		  double *x, double *xex, double *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  double *fwk, int *iparm, double *fparm,
		  int& ier ) {
    dbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}
void pcg::xbasr ( int& ijob, int& ireq,
		  complex<float> *x, complex<float> *xex, complex<float> *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  complex<float> *fwk, int *iparm, complex<float> *fparm,
		  int& ier ) {
    cbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}
void pcg::xbasr ( int& ijob, int& ireq,
		  complex<double> *x, complex<double> *xex, complex<double> *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  complex<double> *fwk, int *iparm, complex<double> *fparm,
		  int& ier ) {
    zbasr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	   ier );
}

//---------------------------------------------------------------------------//
// xgmrsr
//---------------------------------------------------------------------------//

void pcg::xgmrsr( int& ijob, int& ireq,
		  float *x, float *xex, float *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  float *fwk, int *iparm, float *fparm,
		  int& ier ) {
    sgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}
void pcg::xgmrsr( int& ijob, int& ireq,
		  double *x, double *xex, double *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  double *fwk, int *iparm, double *fparm,
		  int& ier ) {
    dgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}
void pcg::xgmrsr( int& ijob, int& ireq,
		  complex<float> *x, complex<float> *xex, complex<float> *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  complex<float> *fwk, int *iparm, complex<float> *fparm,
		  int& ier ) {
    cgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}
void pcg::xgmrsr( int& ijob, int& ireq,
		  complex<double> *x, complex<double> *xex, complex<double> *b,
		  int& iva, int& ivql, int& ivqr, int *iwk,
		  complex<double> *fwk, int *iparm, complex<double> *fparm,
		  int& ier ) {
    zgmrsr( ijob, ireq, x, xex, b, iva, ivql, ivqr, iwk, fwk, iparm, fparm,
	    ier );
}

//---------------------------------------------------------------------------//
//                              end of PCG_Subroutines.cc
//---------------------------------------------------------------------------//
