//----------------------------------*-C++-*----------------------------------//
// PCG_Subroutines.hh
// Dave Nystrom
// Fri May  2 11:02:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_PCG_Subroutines_hh__
#define __linalg_PCG_Subroutines_hh__

#include <complex.h>

#define PCG_NAMESPACE_BEG namespace pcg {
#define PCG_NAMESPACE_END }

PCG_NAMESPACE_BEG

void xdfalt( int *iparm, float           *fparm );
void xdfalt( int *iparm, double          *fparm );
void xdfalt( int *iparm, complex<float>  *fparm );
void xdfalt( int *iparm, complex<double> *fparm );

// Basic iterative method.
void xbasr ( int& ijob, int& ireq,
	     float *x, float *xex, const float *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     float *fwk, int *iparm, float *fparm,
	     int& ier );
void xbasr ( int& ijob, int& ireq,
	     double *x, double *xex, const double *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     double *fwk, int *iparm, double *fparm,
	     int& ier );
void xbasr ( int& ijob, int& ireq,
	     complex<float> *x, complex<float> *xex, const complex<float> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<float> *fwk, int *iparm, complex<float> *fparm,
	     int& ier );
void xbasr ( int& ijob, int& ireq,
	     complex<double> *x, complex<double> *xex, const complex<double> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<double> *fwk, int *iparm, complex<double> *fparm,
	     int& ier );

// Restarted GMRES.
void xgmrsr( int& ijob, int& ireq,
	     float *x, float *xex, const float *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     float *fwk, int *iparm, float *fparm,
	     int& ier );
void xgmrsr( int& ijob, int& ireq,
	     double *x, double *xex, const double *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     double *fwk, int *iparm, double *fparm,
	     int& ier );
void xgmrsr( int& ijob, int& ireq,
	     complex<float> *x, complex<float> *xex, const complex<float> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<float> *fwk, int *iparm, complex<float> *fparm,
	     int& ier );
void xgmrsr( int& ijob, int& ireq,
	     complex<double> *x, complex<double> *xex, const complex<double> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<double> *fwk, int *iparm, complex<double> *fparm,
	     int& ier );

// Conjugate gradient.
void xcgr  ( int& ijob, int& ireq,
	     float *x, float *xex, const float *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     float *fwk, int *iparm, float *fparm,
	     int& ier );
void xcgr  ( int& ijob, int& ireq,
	     double *x, double *xex, const double *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     double *fwk, int *iparm, double *fparm,
	     int& ier );
void xcgr  ( int& ijob, int& ireq,
	     complex<float> *x, complex<float> *xex, const complex<float> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<float> *fwk, int *iparm, complex<float> *fparm,
	     int& ier );
void xcgr  ( int& ijob, int& ireq,
	     complex<double> *x, complex<double> *xex, const complex<double> *b,
	     int& iva, int& ivql, int& ivqr, int *iwk,
	     complex<double> *fwk, int *iparm, complex<double> *fparm,
	     int& ier );

// Enumerations.

enum { BASIC = 1, BCGSTAB, BCGSTAB2, CGS, TFQMR, GMRES, GMRES_H, OMIN,
       ORES, IOM, CG, BCG, CGNE, CGNP, CGNR, BCGSTABL, GMRES_CH, GMRES_R,
       LRES, LSQE, LSQR, RBCG, SQMR };
enum { NOUT = 1, LEVOUT, NRU, ITSMAX, ITS, MALLOC, NWI, NWF, NWIUSD,
       NWFUSD, IPTR, NTEST, IQSIDE, IUINIT, NEEDRC, NS1, NS2, ICKSTG,
       IUEXAC, IDOT, ISTATS, ITIMER, ICOMM, MSGMIN, MSGMAX, MSGTYP,
       ICLEV, ISCALE = 41, ICTRAN, ILDIM };
enum { CTIMER = 1, RTIMER, FLOPSR, ZETA, STPTST, ALPHA, RELRSD, RELERR,
       CTIMEI = 11, RTIMEI, FLOPSI, CTIMET = 21, RTIMET, FLOPST, OMEGA };
enum { WRNUNK = 1, WRNNOC, WRNZET, WRNSTT };
enum { ERRMVP = -17, ERRPRC, ERRARR, ERRMEM, ERRIO, ERRCOM, ERRMAT,
       ERRSTP, ERRMDI, ERRSTG, ERRBDB, ERRBDI, ERRINA, ERRARG, ERRIFW,
       ERRIIW, ERRUNK };
enum { NIPARM = 50, NFPARM = 30 };
enum { YES = 1, NO = 0, DFALT = -1 };
enum { QNONE, QLEFT, QRIGHT, QSPLIT };
enum { USZERO = -2, UZERO = 0, UNZERO, USRAND, UPRAND };
enum { LEV0, LEVERR, LEVWRN, LEVIT, LEVPRM, LEVALG };
enum { TSTUSR = -3, TSTEX, TSTDFA, TST0, TSTSE, TSTSR, TSTSLR, TSTSRR,
       TSTRE, TSTRR, TSTRLR, TSTRRR };
enum { JTERM = -1, JIRT, JINIT, JINITA, JRUN, JRUNA, JRUNAQ, JTEST = 9 };
enum { JAV = 3, JATV, JQLV, JQLTV, JQRV, JQRTV };
enum { SCNONE, SCSD, SCLD, SCLR };
enum { GRBIN, GRGRAY };
enum { ORRG, ORGR };
enum { ORRCSG, ORRCGS, ORRSCG, ORRSGC, ORRGCS, ORRGSC, ORCRSG, ORCRGS,
       ORCSRG, ORCSGR, ORCGRS, ORCGSR, ORSRCG, ORSRGC, ORSCRG, ORSCGR,
       ORSGRC, ORSGCR, ORGRCS, ORGRSC, ORGCRS, ORGCSR, ORGSRC, ORGSCR };
enum { FALSE, TRUE };

PCG_NAMESPACE_END

#undef PCG_NAMESPACE_BEG
#undef PCG_NAMESPACE_END

// Handle machine specific Fortran name linkage.

#if defined(sun) || defined(__sgi) || defined(__linux)
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

extern "C" {
    void sdfalt( int *iparm, float           *fparm );
    void ddfalt( int *iparm, double          *fparm );
    void cdfalt( int *iparm, complex<float>  *fparm );
    void zdfalt( int *iparm, complex<double> *fparm );

// Basic iterative method.
    void sbasr ( int& ijob, int& ireq,
		 float *x, float *xex, const float *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 float *fwk, int *iparm, float *fparm,
		 int& ier );
    void dbasr ( int& ijob, int& ireq,
		 double *x, double *xex, const double *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 double *fwk, int *iparm, double *fparm,
		 int& ier );
    void cbasr ( int& ijob, int& ireq,
		 complex<float> *x, complex<float> *xex, const complex<float> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<float> *fwk, int *iparm, complex<float> *fparm,
		 int& ier );
    void zbasr ( int& ijob, int& ireq,
		 complex<double> *x, complex<double> *xex, const complex<double> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<double> *fwk, int *iparm, complex<double> *fparm,
		 int& ier );

// Restarted GMRES.
    void sgmrsr( int& ijob, int& ireq,
		 float *x, float *xex, const float *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 float *fwk, int *iparm, float *fparm,
		 int& ier );
    void dgmrsr( int& ijob, int& ireq,
		 double *x, double *xex, const double *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 double *fwk, int *iparm, double *fparm,
		 int& ier );
    void cgmrsr( int& ijob, int& ireq,
		 complex<float> *x, complex<float> *xex, const complex<float> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<float> *fwk, int *iparm, complex<float> *fparm,
		 int& ier );
    void zgmrsr( int& ijob, int& ireq,
		 complex<double> *x, complex<double> *xex, const complex<double> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<double> *fwk, int *iparm, complex<double> *fparm,
		 int& ier );

// Conjugate gradient.
    void scgr  ( int& ijob, int& ireq,
		 float *x, float *xex, const float *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 float *fwk, int *iparm, float *fparm,
		 int& ier );
    void dcgr  ( int& ijob, int& ireq,
		 double *x, double *xex, const double *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 double *fwk, int *iparm, double *fparm,
		 int& ier );
    void ccgr  ( int& ijob, int& ireq,
		 complex<float> *x, complex<float> *xex, const complex<float> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<float> *fwk, int *iparm, complex<float> *fparm,
		 int& ier );
    void zcgr  ( int& ijob, int& ireq,
		 complex<double> *x, complex<double> *xex, const complex<double> *b,
		 int& iva, int& ivql, int& ivqr, int *iwk,
		 complex<double> *fwk, int *iparm, complex<double> *fparm,
		 int& ier );
}

#endif                          // __linalg_PCG_Subroutines_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_Subroutines.hh
//---------------------------------------------------------------------------//
