//----------------------------------*-C++-*----------------------------------//
// PCG_Ctrl.hh
// Dave Nystrom
// Mon Jan 13 17:40:28 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_PCG_Ctrl_hh__
#define __linalg_PCG_Ctrl_hh__

#include "ds++/Mat.hh"
#include "ds++/SP.hh"

#include "PCG_MatVec.hh"
#include "PCG_PreCond.hh"

//===========================================================================//
// class PCG_Ctrl - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original by Dave Nystrom
// 1) extensively re-vamped by Rob Lowrie
// 
//===========================================================================//

template<class T>
class PCG_Ctrl
{
  public:

    // PUBLIC TYPES

    // Iterative methods (not all PCG methods are supported)
    enum Method { BAS, CG, GMRS };

    // Indices in IPARMS array
    enum Iparms { NOUT = 1, LEVOUT, NRU, ITSMAX, ITS, MALLOC, NWI, NWF, NWIUSD,
		  NWFUSD, IPTR, NTEST, IQSIDE, IUINIT, NEEDRC, NS1, NS2, ICKSTG,
		  IUEXAC, IDOT, ISTATS, ITIMER, ICOMM, MSGMIN, MSGMAX, MSGTYP,
		  ISCALE = 41, ICTRAN };

    // Indices in FPARMS array
    enum Fparms { CTIMER = 1, RTIMER, FLOPSR, ZETA, STPTST, ALPHA, RELRSD,
		  RELERR, CTIMEI = 11, RTIMEI, FLOPSI, CTIMET = 21, RTIMET,
		  FLOPST, OMEGA };

    // Legal values for LEVOUT in IPARMS (output level)
    enum OutputLevel { LEV0, LEVERR, LEVWRN, LEVIT, LEVPRM, LEVALG };

    // Legal values for NTEST in IPARMS (stop tests)
    enum StopTest { TSTUSR = -3, TSTEX, TSTDFA, TST0, TSTSE, TSTSR, TSTSLR,
		    TSTSRR, TSTRE, TSTRR, TSTRLR, TSTRRR };

    // Legal values for NTEST in IPARMS (preconditioners)
    enum Precon { QNONE, QLEFT, QRIGHT, QSPLIT };

    // Legal values for IUINIT in IPARMS (how initial guess is set)
    enum Uinit { USZERO = -2, UDFALT = -1, UZERO = 0, UNZERO, USRAND, UPRAND };

    // Logic values
    enum Logical { DFALT = -1, NO = 0, YES = 1 };

  private:

    // PRIVATE TYPES

    // Warning codes (not used yet)
    //    enum Warnings { WRNUNK = 1, WRNNOC, WRNZET, WRNSTT };

    // IJOB values
    enum { JTERM = -1, JIRT, JINIT, JINITA, JRUN, JRUNA, JRUNAQ, JTEST = 9 };

    // IREQ values
    enum { JAV = 3, JATV, JQLV, JQLTV, JQRV, JQRTV };

    // DATA

    dsxx::Mat1<int> d_iparm; // PCG IPARM array
    dsxx::Mat1<T> d_fparm; // PCG FPARM array
    dsxx::Mat1<T> d_uExact; // PCG UEXACT array
    Method d_method; // The PCG iterative method to be called

  public:

    // CREATORS
    
    PCG_Ctrl(const Method method = GMRS);
    PCG_Ctrl(const PCG_Ctrl &rhs);
    ~PCG_Ctrl() {}

    // ACCESSORS

    int getIparm(const Iparms parm) const { return d_iparm(parm); }
    T getFparm(const Fparms parm) const { return d_fparm(parm); }
    void printParams() const;

    // MANIPULATORS

    PCG_Ctrl& operator=(const PCG_Ctrl &rhs);
    void solve(dsxx::Mat1<T>& x,
	       const dsxx::Mat1<T>& b,
	       dsxx::SP< PCG_MatVec<T> > pcg_matvec,
	       dsxx::SP< PCG_PreCond<T> > pcg_precond,
	       const int nru = -1);
    void setIparm(const Iparms parm,
		  const int value);
    void setFparm(const Fparms parm,
		  const T value);
    void setLogical(const Iparms parm,
		    const Logical value);
    void setOutputLevel(const OutputLevel value);
    void setStopTest(const StopTest value);
    void setPrecon(const Precon value);
    void setUinit(const Uinit value);
    void setUexact(const dsxx::Mat1<T>& uExact);

  private:

    // IMPLEMENTATION

    void computeWorkSpace();
    void callPCG(dsxx::Mat1<T>& x,
		 const dsxx::Mat1<T>& b,
		 int &ijob,
		 int &ireq,
		 int &iva,
		 int &ivql,
		 int &ivqr,
		 dsxx::Mat1<int> &iwork,
		 dsxx::Mat1<T> &fwork);
    int getSize() const { return d_iparm(NRU); }
};

#endif                          // __linalg_PCG_Ctrl_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_Ctrl.hh
//---------------------------------------------------------------------------//
