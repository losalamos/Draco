//----------------------------------*-C++-*----------------------------------//
// PCG_Ctrl.cc
// Dave Nystrom
// Mon Jan 13 17:40:29 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "PCG_Ctrl.hh"
#include "PCG_Subroutines.hh"
#include <iostream>
#include <sstream>

namespace rtt_pcgWrap {

//---------------------------------------------------------------------------//
// Constructor.
//---------------------------------------------------------------------------//

template<class T>
PCG_Ctrl<T>::
PCG_Ctrl(const Method method)
    : d_iparm(rtt_dsxx::Bounds(1,50)),
      d_fparm(rtt_dsxx::Bounds(1,30)),
      d_uExact(1),
      d_method(method)
{
    // is_supported must be checked for all constructors.
    Insist(is_supported(), "PCG unsupported on this platform!");

    // Initialize iparm and fparm arrays via PCG defaults
    xdfalt(&d_iparm(1), &d_fparm(1));
    
    d_iparm(MALLOC) = NO;  // Don't allow this for now; see allocateWorkArrays
}

//---------------------------------------------------------------------------//
// Copy constructor.
//---------------------------------------------------------------------------//

template<class T>
PCG_Ctrl<T>::
PCG_Ctrl(const PCG_Ctrl<T> &rhs)
    : d_iparm(rhs.d_iparm),
      d_fparm(rhs.d_fparm),
      d_uExact(rhs.d_uExact),
      d_method(rhs.d_method)
{
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

template<class T>
PCG_Ctrl<T> &
PCG_Ctrl<T>::
operator=(const PCG_Ctrl<T> &rhs)
{
    if ( this == &rhs ) {
	return *this;
    }

    d_iparm = rhs.d_iparm;
    d_fparm = rhs.d_fparm;
    d_uExact = rhs.d_uExact;
    d_method = rhs.d_method;

    return *this;
}

//---------------------------------------------------------------------------//
// Allocates the required workspace.
//
// nru: Length of solution vector.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
allocateWorkspace(const int nru)
{
    Require(nru > 0);

    d_iparm(NRU) = nru;
    
    computeWorkSpace();

    d_iwork.redim(d_iparm(NWI));
    d_fwork.redim(d_iparm(NWF));
}

//---------------------------------------------------------------------------//
// Main controller method.
//
// x: Output solution and possibly input initial guess
// b: right-hand size vector
// pcg_matvec: matvec routine
// pcg_precond: preconditioner routine
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
solve(rtt_dsxx::Mat1<T>& x,
      const rtt_dsxx::Mat1<T>& b,
      rtt_dsxx::SP< PCG_MatVec<T> > pcg_matvec,
      rtt_dsxx::SP< PCG_PreCond<T> > pcg_precond)
{
    if ( d_iwork.size() == 0 && d_fwork.size() == 0 ) {
	// allocateWorkspace has not been called

	Ensure(x.size() == b.size()); // vague otherwise
	allocateWorkspace(x.size());
    }
    else {
	Ensure(d_iparm(NRU) <= static_cast<int>(x.size()));
	Ensure(d_iparm(NRU) <= static_cast<int>(b.size()));
    }

    using rtt_dsxx::Mat1;

    // These variables are used to communicate with the PCG package.
    // See the "call _methR" PCG documentation.
    
    int ijob;
    int ireq;
    int iva;
    int ivql;
    int ivqr;
    
    // Initialize ijob so that initialization is done on the first
    // callPCG call.
    
    ijob = JINIT;

    // Loop on PCG calls in the reverse communication mode.
    
    bool done = false;

    while ( ! done ) {

	callPCG(x, b, ijob, ireq, iva, ivql, ivqr);

	ijob = JRUN;  // from now on, RUN!

	// Handle the request from callPCG

	switch ( ireq ) {

	case JTERM: {
	    done = true;
	    break;
	}

	case JAV: {
	    // cout << "Preparing for MatVec." << endl << flush;
	    const Mat1<T> xmatvec(&d_fwork(ivqr-1), getSize());
	    Mat1<T> bmatvec(&d_fwork(iva-1), getSize());
	    pcg_matvec->MatVec(bmatvec, xmatvec);
	    // cout << "Done with     MatVec." << endl << flush;
	    break;
	}
	
	case JQLV: {
	    // cout << "Preparing for Left_PreCond." << endl << flush;
	    Mat1<T> xprecond(&d_fwork(ivql-1), getSize());
	    const Mat1<T> bprecond(&d_fwork(iva-1), getSize());
	    pcg_precond->Left_PreCond(xprecond, bprecond);
	    // cout << "Done with     Left_PreCond." << endl << flush;
	    break;
	}
	
	case JQRV: {
	    // cout << "Preparing for Right_PreCond." << endl << flush;
	    Mat1<T> xprecond(&d_fwork(ivqr-1), getSize());
	    const Mat1<T> bprecond(&d_fwork(ivql-1), getSize());
	    pcg_precond->Right_PreCond(xprecond, bprecond);
	    // cout << "Done with     Right_PreCond." << endl << flush;
	    break;
	}

	default: {
	    std::ostringstream mesg;
	    mesg << "PCG returned IREQ = " << ireq
		 << " which PCG_Ctrl::solve cannot handle.";
	    Insist(0, mesg.str().c_str());
	}
	}
    } // end of while ( ! done )
}

//---------------------------------------------------------------------------//
// Call a pcg iterative method in reverse communication mode.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
callPCG(rtt_dsxx::Mat1<T> &x,
	const rtt_dsxx::Mat1<T> &b,
	int &ijob,
	int &ireq,
	int &iva,
	int &ivql,
	int &ivqr)
{
    int iError = 0;
    
    switch ( d_method ) {
    case BAS: {
	xbasr(ijob, ireq, &x(0), &d_uExact(0), &b(0), iva, ivql,
	      ivqr, &d_iwork(0), &d_fwork(0), &d_iparm(1), &d_fparm(1),
	      iError);
	break;
    }
    
    case CG: {
	xcgr(ijob, ireq, &x(0), &d_uExact(0), &b(0), iva, ivql,
	     ivqr, &d_iwork(0), &d_fwork(0), &d_iparm(1), &d_fparm(1),
	     iError);
	break;
    }
    
    case GMRS: {
	xgmrsr(ijob, ireq, &x(0), &d_uExact(0), &b(0), iva, ivql,
	       ivqr, &d_iwork(0), &d_fwork(0), &d_iparm(1), &d_fparm(1),
	       iError);
	break;
    }
    }

    Ensure(iError == 0);
}

//---------------------------------------------------------------------------//
// This routine sets the memory parameters NWI and NWF in the IPARM array,
// depending on the various methods.  This is a translation of the PCG
// Reference Manual, Chapter 9.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
computeWorkSpace()
{
    int nwi;
    int nwf;
    
    if( d_iparm(MALLOC) == YES ) {
	// PCG mallocs its own memory
	// For the present constructor, we don't allow this, but leave
	// this logic here for now.
	nwi = 1;
	nwf = 1;
    }
    else {
	// PCG actually uses the work arrays

	// First translate Table 9.1 of the PCG RM.  Ignore CMF case.
	// Since we're never actually passing PCG a matrix, we can
	// also ignore the GR format stuff.  Therefore, SPMD and Uniprocessor
	// are the same (phew!)
	const int nwiv = 0;
	const int nwfv = d_iparm(NRU) + 2;
	const int nwigenl = 32 + 2 * nwiv;
	const int nwfgenl = 32 + 2 * nwfv;
	
	int nwistat;
	if ( d_iparm(ISTATS) == NO ) {
	    nwistat = 0;
	}
	else if ( d_iparm(ISTATS) == YES ) {
	    nwistat = 20 + 4 * nwfv;
	}
	else {
	    Insist(0, "ISTATS has illegal value.");
	}

	int nwitst;
	if ( d_iparm(NTEST) == TST0 ||
	     d_iparm(NTEST) == TSTDFA ) {
	    nwitst = 0;
	}
	else {
	    nwitst = 2 * nwfv;
	}

	// Translate Table 9.2

	int nwiit;
	int nwfit;
	// const int ns1 = d_iparm(NS1);
	const int ns2 = d_iparm(NS2);

	switch ( d_method ) {
	case BAS:
	    nwiit = 26 + 5 * nwiv;
	    nwfit = 7  + 5 * nwfv;
	    break;

	case CG:
	    nwiit = 25 + 5 * nwiv;
	    nwfit = 12 + 5 * nwfv;
	    break;

	case GMRS: {
	    // Assume that NR = 0 for now
	    if ( d_iparm(IQSIDE) == QNONE ||
		 d_iparm(IQSIDE) == QLEFT ) {
		nwiit = 39 + (ns2 + 4) * nwiv;
		nwfit = ns2 * (ns2 + 9) + 31 + (ns2 + 4) * nwfv;
	    }
	    else {
		nwiit = 39 + (3 * ns2 + 8) * nwiv;
		nwfit = ns2 * (ns2 + 9) + 31 + (3 * ns2 + 8) * nwfv;
	    }

	    break;
	}
	}

	// Sum up the requirements.  Again, nwiscl = nwfscl = 0 since we
	// don't use GR format.  nwiprec = nwfprec = nwfstat = nwftst = 0,
	// since we have no idea what these are.

	nwi = nwigenl + nwiv + nwiit + nwistat + nwitst;
	nwf = nwfgenl + nwfv + nwfit;
    }

    d_iparm(NWI) = nwi;
    d_iparm(NWF) = nwf;
}

//---------------------------------------------------------------------------//
// Set a value for IPARM array.  Only certain values are permitted to be
// changed.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setIparm(const Iparms parm,
	 const int value)
{
    switch ( parm ) {
    case NOUT:
    case ITSMAX:
    case NS1:
    case NS2:
    case ITIMER: // Do we really need Connection Machine stuff???
    case ICOMM:
    case MSGMIN:
    case MSGMAX:
    case MSGTYP: {
	d_iparm(parm) = value;
	break;
    }
    default: {
	std::ostringstream mesg;
	mesg << "Cannot use PCG_Ctrl::setIparm to set IPARM index = "
	     << parm;
	Insist(0, mesg.str().c_str());
	break;
    }
    }
}

//---------------------------------------------------------------------------//
// Set a value for FPARM array.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setFparm(const Fparms parm,
	 const T value)
{
    d_fparm(parm) = value;
}

//---------------------------------------------------------------------------//
// Set the exact solution.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setUexact(const rtt_dsxx::Mat1<T>& uExact)
{
    d_uExact = uExact;
    d_iparm(IUEXAC) = YES;
}

//---------------------------------------------------------------------------//
// Set a logical value in IPARM array.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setLogical(const Iparms parm,
	   const Logical value)
{
    switch ( parm ) {
    case ICKSTG:
    case IDOT:
    case ISTATS: {
	d_iparm(parm) = value;
	break;
    }
    default: {
	std::ostringstream mesg;
	mesg << "Cannot use PCG_Ctrl::setLogical to set IPARM index = "
	     << parm;
	Insist(0, mesg.str().c_str());
	break;
    }
    }
}

//---------------------------------------------------------------------------//
// Set the output level in the IPARM array
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setOutputLevel(const OutputLevel value)
{
    d_iparm(LEVOUT) = value;
}

//---------------------------------------------------------------------------//
// Set the stopping test in the IPARM array
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setStopTest(const StopTest value)
{
    d_iparm(NTEST) = value;
}

//---------------------------------------------------------------------------//
// Set the preconditioner in the IPARM array
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setPrecon(const Precon value)
{
    d_iparm(IQSIDE) = value;
}

//---------------------------------------------------------------------------//
// Set how initial condition is set in the IPARM array
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
setUinit(const Uinit value)
{
    d_iparm(IUINIT) = value;
}

//---------------------------------------------------------------------------//
// Print values for PCG iparm and fparm arrays.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::
printParams() const
{
    using std::cout;
    using std::endl;
    
// Revcom level parameters.
    cout << "----------------------------------------------" << endl;
    cout << "Revcom level parameters."                       << endl;
    cout << "----------------------------------------------" << endl;
    cout << "     nout   = " << d_iparm(NOUT)   << endl;
    cout << "     levout = " << d_iparm(LEVOUT) << endl;
    cout << "     nru    = " << d_iparm(NRU)    << endl;
    cout << "     itsmax = " << d_iparm(ITSMAX) << endl;
    cout << "     its    = " << d_iparm(ITS)    << endl;
    cout << "     malloc = " << d_iparm(MALLOC) << endl;
    cout << "     nwi    = " << d_iparm(NWI)    << endl;
    cout << "     nwf    = " << d_iparm(NWF)    << endl;
    cout << "     nwiusd = " << d_iparm(NWIUSD) << endl;
    cout << "     nwfusd = " << d_iparm(NWFUSD) << endl;
    cout << "     iptr   = " << d_iparm(IPTR)   << endl;
    cout << "     ntest  = " << d_iparm(NTEST)  << endl;
    cout << "     iqside = " << d_iparm(IQSIDE) << endl;
    cout << "     iuinit = " << d_iparm(IUINIT) << endl;
    cout << "     needrc = " << d_iparm(NEEDRC) << endl;
    cout << "     ns1    = " << d_iparm(NS1)    << endl;
    cout << "     ns2    = " << d_iparm(NS2)    << endl;
    cout << "     ickstg = " << d_iparm(ICKSTG) << endl;
    cout << "     iuexac = " << d_iparm(IUEXAC) << endl;
    cout << "     idot   = " << d_iparm(IDOT)   << endl;
    cout << "     istats = " << d_iparm(ISTATS) << endl;
    cout << "     itimer = " << d_iparm(ITIMER) << endl;
    cout << "     icomm  = " << d_iparm(ICOMM)  << endl;
    cout << "     msgmin = " << d_iparm(MSGMIN) << endl;
    cout << "     msgmax = " << d_iparm(MSGMAX) << endl;
    cout << "     msgtyp = " << d_iparm(MSGTYP) << endl;
    cout << " "                                 << endl;
    cout << "     ctimer = " << d_fparm(CTIMER) << endl;
    cout << "     rtimer = " << d_fparm(RTIMER) << endl;
    cout << "     flopsr = " << d_fparm(FLOPSR) << endl;
    cout << "     zeta   = " << d_fparm(ZETA)   << endl;
    cout << "     stptst = " << d_fparm(STPTST) << endl;
    cout << "     alpha  = " << d_fparm(ALPHA)  << endl;
    cout << "     relrsd = " << d_fparm(RELRSD) << endl;
    cout << "     relerr = " << d_fparm(RELERR) << endl;
}

} // namespace rtt_pcgWrap

//---------------------------------------------------------------------------//
//                              end of PCG_Ctrl.cc
//---------------------------------------------------------------------------//
