//----------------------------------*-C++-*----------------------------------//
// PCG_Ctrl.cc
// Dave Nystrom
// Mon Jan 13 17:40:29 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/PCG_Ctrl.hh"
#include "linalg/PCG_Subroutines.hh"

//---------------------------------------------------------------------------//
// Constructor.
//---------------------------------------------------------------------------//

template<class T>
PCG_Ctrl<T>::PCG_Ctrl( const pcg_DB& pcg_db, int _nru )
    : pcg_DB(pcg_db),
      iparm(Bounds(1,50)), fparm(Bounds(1,30)),
      iwk(pcg_db.nwi), fwk(pcg_db.nwf),
      xex(1), nru(_nru)
{
// Initialize some stuff from pcg_DB.
    itmeth = pcg_db.itmeth;

// Compute required size of iwk and fwk and resize them.
    set_nwi();
    set_nwf();

    iwk.redim( nwi );
    fwk.redim( nwf );

// Initialize iparm and fparm arrays.
    set_default();

    iparm(pcg::NOUT)   = nout;
    iparm(pcg::LEVOUT) = levout;
    iparm(pcg::NRU)    = nru;
    iparm(pcg::ITSMAX) = itsmax;
    iparm(pcg::MALLOC) = malloc;
    iparm(pcg::NWI)    = nwi;
    iparm(pcg::NWF)    = nwf;
    iparm(pcg::NTEST)  = ntest;
    iparm(pcg::IQSIDE) = iqside;
    iparm(pcg::IUINIT) = iuinit;
    iparm(pcg::NEEDRC) = needrc;
    iparm(pcg::NS1)    = ns1;
    iparm(pcg::NS2)    = ns2;
    iparm(pcg::ICKSTG) = ickstg;
    iparm(pcg::IUEXAC) = iuexac;
    iparm(pcg::IDOT)   = idot;
    iparm(pcg::ISTATS) = istats;

    fparm(pcg::CTIMER) = ctimer;
    fparm(pcg::RTIMER) = rtimer;
    fparm(pcg::FLOPSR) = flopsr;
    fparm(pcg::ZETA)   = zeta;
    fparm(pcg::ALPHA)  = alpha;

// Print control parameters.
//    print_params();
}

//---------------------------------------------------------------------------//
// Main controller method.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::pcg_fe( Mat1<T>& x, Mat1<T>& b,
			  SP< PCG_MatVec<T> >& pcg_matvec,
			  SP< PCG_PreCond<T> >& pcg_precond )
{
// Initialize ijob.
    ijob = pcg::JINIT;

// Call an iterative method.
    int done=0;
    while (!done) {
	it_method( x, b, xex );

	ijob = pcg::JRUN;

	if( ireq == pcg::JTERM ) {
	    done = 1;
	}
	else if( ireq == pcg::JAV ) {
	//	    cout << "Preparing for MatVec." << endl << flush;
	    Mat1<T> xmatvec(&fwk(ivqr-1),nru);
	    Mat1<T> bmatvec(&fwk(iva-1), nru);
	    pcg_matvec->MatVec(bmatvec,xmatvec);
	//	    cout << "Done with     MatVec." << endl << flush;
	}
	else if( ireq == pcg::JQLV ) {
	//	    cout << "Preparing for Left_PreCond." << endl << flush;
	    Mat1<T> xprecond(&fwk(ivql-1),nru);
	    Mat1<T> bprecond(&fwk(iva-1), nru);
	    pcg_precond->Left_PreCond(xprecond,bprecond);
	//	    cout << "Done with     Left_PreCond." << endl << flush;
	}
	else if( ireq == pcg::JQRV ) {
	    Mat1<T> xprecond(&fwk(ivqr-1),nru);
	    Mat1<T> bprecond(&fwk(ivql-1),nru);
	    pcg_precond->Right_PreCond(xprecond,bprecond);
	}
	else if( ireq == pcg::JTEST ) {
	}
	else if( ireq == pcg::JATV ) {
	}
	else if( ireq == pcg::JQLTV ) {
	}
	else if( ireq == pcg::JQRTV ) {
	}
    }
}

//---------------------------------------------------------------------------//
// Set default values for PCG iparm and fparm arrays.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::set_default()
{
    pcg::xdfalt( &iparm(1), &fparm(1) );
}

//---------------------------------------------------------------------------//
// Call a pcg iterative method.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::it_method( Mat1<T>& x, Mat1<T>& b, Mat1<T>& xex )
{
    if( itmeth == pcg::BASIC ) {
	pcg::xbasr( ijob, ireq, &x(0), &xex(0), &b(0), iva, ivql, ivqr,
		    &iwk(0), &fwk(0), &iparm(1), &fparm(1), ier );
	imatvec = pcg::TRUE;
    }
    else if( itmeth == pcg::GMRES ) {
	pcg::xgmrsr( ijob, ireq, &x(0), &xex(0), &b(0), iva, ivql, ivqr,
		     &iwk(0), &fwk(0), &iparm(1), &fparm(1), ier );
	imatvec = pcg::TRUE;
    }
    else if( itmeth == pcg::CG ) {
	pcg::xcgr( ijob, ireq, &x(0), &xex(0), &b(0), iva, ivql, ivqr,
		   &iwk(0), &fwk(0), &iparm(1), &fparm(1), ier );
	imatvec = pcg::TRUE;
    }
    else {
	throw("Need to choose a valid pcg iterative method.");
    }
}

//---------------------------------------------------------------------------//
// Set nwi.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::set_nwi()
{
    if( itmeth == pcg::BASIC    ) nwi = 100;
    if( itmeth == pcg::BCGSTAB  ) nwi = 100;
    if( itmeth == pcg::BCGSTAB2 ) nwi = 100;
    if( itmeth == pcg::BCGSTABL ) nwi = 100;
    if( itmeth == pcg::CGS      ) nwi = 100;
    if( itmeth == pcg::TFQMR    ) nwi = 100;
    if( itmeth == pcg::GMRES    ) nwi = 100;
    if( itmeth == pcg::GMRES_H  ) nwi = 100;
    if( itmeth == pcg::OMIN     ) nwi = 100;
    if( itmeth == pcg::ORES     ) nwi = 100;
    if( itmeth == pcg::IOM      ) nwi = 100;
    if( itmeth == pcg::CG       ) nwi = 100;
    if( itmeth == pcg::BCG      ) nwi = 100;

    if( malloc == pcg::TRUE     ) nwi = 1;
}

//---------------------------------------------------------------------------//
// Set nwf.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::set_nwf()
{
    int nrup2   = nru + 2;
    int nwfgenl = 31 + nrup2*2;
    int nwfit   = 0;
    int nwfstat = 0;
    int nwftst  = 0;

    if( itmeth == pcg::BASIC    ) nwfit =  7 + nrup2*5;
    if( itmeth == pcg::BCGSTAB  ) nwfit = 15 + nrup2*13;
    if( itmeth == pcg::BCGSTAB2 ) nwfit = 15 + nrup2*25;
    if( itmeth == pcg::BCGSTABL ) nwfit = 31 + nrup2*(2*ns2+8) + ns2*9
				        + ns2*ns2;
    if( itmeth == pcg::CGS      ) nwfit = 12 + nrup2*11;
    if( itmeth == pcg::TFQMR    ) nwfit = 16 + nrup2*18;
    if( itmeth == pcg::CG       ) nwfit = 12 + nrup2*5;
    if( itmeth == pcg::BCG      ) nwfit = 12 + nrup2*9;
    if( iqside >= pcg::QRIGHT ) {
	if( itmeth == pcg::GMRES   ) nwfit = 31 + nrup2*(2*ns2+8) + ns2*9
					   + ns2*ns2;
	if( itmeth == pcg::GMRES_H ) nwfit = 26 + nrup2*(2*ns2+7) + ns2*7
					   + ns2*ns2;
	if( itmeth == pcg::OMIN    ) nwfit = 13 + nrup2*(3*ns1+5) + ns1;
	if( itmeth == pcg::ORES    ) nwfit = 12 + nrup2*(4*ns1+6) + ns1;
	if( itmeth == pcg::IOM     ) nwfit = 31 + nrup2*(3*ns1+9) + ns1*5;
    }        
    else {
	if( itmeth == pcg::GMRES   ) nwfit = 31 + nrup2*(1*ns2+4) + ns2*9
					   + ns2*ns2;
	if( itmeth == pcg::GMRES_H ) nwfit = 26 + nrup2*(2*ns2+7) + ns2*7
					   + ns2*ns2;
	if( itmeth == pcg::OMIN    ) nwfit = 13 + nrup2*(2*ns1+3) + ns1;
	if( itmeth == pcg::ORES    ) nwfit = 12 + nrup2*(2*ns1+4) + ns1;
	if( itmeth == pcg::IOM     ) nwfit = 31 + nrup2*(2*ns1+7) + ns1*5;
    }

    if( istats == 1 ) nwfstat = 20 + nrup2*4;

    if( ntest != pcg::TST0 && ntest != pcg::TSTDFA ) {
	nwftst = nrup2*2;
    }

    if( malloc == pcg::TRUE ) {
	nwf = 1;
    }
    else {
	nwf = nwfgenl + nwfstat + nwftst + nwfit;
    }
}

//---------------------------------------------------------------------------//
// Print values for PCG iparm and fparm arrays.
//---------------------------------------------------------------------------//

template<class T>
void PCG_Ctrl<T>::print_params()
{
// Revcom level parameters.
    cout << "----------------------------------------------" << endl;
    cout << "Revcom level parameters."                       << endl;
    cout << "----------------------------------------------" << endl;
    cout << "     nout   = " << iparm(pcg::NOUT)   << endl;
    cout << "     levout = " << iparm(pcg::LEVOUT) << endl;
    cout << "     nru    = " << iparm(pcg::NRU)    << endl;
    cout << "     itsmax = " << iparm(pcg::ITSMAX) << endl;
    cout << "     its    = " << iparm(pcg::ITS)    << endl;
    cout << "     malloc = " << iparm(pcg::MALLOC) << endl;
    cout << "     nwi    = " << iparm(pcg::NWI)    << endl;
    cout << "     nwf    = " << iparm(pcg::NWF)    << endl;
    cout << "     nwiusd = " << iparm(pcg::NWIUSD) << endl;
    cout << "     nwfusd = " << iparm(pcg::NWFUSD) << endl;
    cout << "     iptr   = " << iparm(pcg::IPTR)   << endl;
    cout << "     ntest  = " << iparm(pcg::NTEST)  << endl;
    cout << "     iqside = " << iparm(pcg::IQSIDE) << endl;
    cout << "     iuinit = " << iparm(pcg::IUINIT) << endl;
    cout << "     needrc = " << iparm(pcg::NEEDRC) << endl;
    cout << "     ns1    = " << iparm(pcg::NS1)    << endl;
    cout << "     ns2    = " << iparm(pcg::NS2)    << endl;
    cout << "     ickstg = " << iparm(pcg::ICKSTG) << endl;
    cout << "     iuexac = " << iparm(pcg::IUEXAC) << endl;
    cout << "     idot   = " << iparm(pcg::IDOT)   << endl;
    cout << "     istats = " << iparm(pcg::ISTATS) << endl;
    cout << "     itimer = " << iparm(pcg::ITIMER) << endl;
    cout << "     icomm  = " << iparm(pcg::ICOMM)  << endl;
    cout << "     msgmin = " << iparm(pcg::MSGMIN) << endl;
    cout << "     msgmax = " << iparm(pcg::MSGMAX) << endl;
    cout << "     msgtyp = " << iparm(pcg::MSGTYP) << endl;
    cout << "     iclev  = " << iparm(pcg::ICLEV)  << endl;
    cout << " "                                    << endl;
    cout << "     ctimer = " << fparm(pcg::CTIMER) << endl;
    cout << "     rtimer = " << fparm(pcg::RTIMER) << endl;
    cout << "     flopsr = " << fparm(pcg::FLOPSR) << endl;
    cout << "     zeta   = " << fparm(pcg::ZETA)   << endl;
    cout << "     stptst = " << fparm(pcg::STPTST) << endl;
    cout << "     alpha  = " << fparm(pcg::ALPHA)  << endl;
    cout << "     relrsd = " << fparm(pcg::RELRSD) << endl;
    cout << "     relerr = " << fparm(pcg::RELERR) << endl;
}

//---------------------------------------------------------------------------//
//                              end of PCG_Ctrl.cc
//---------------------------------------------------------------------------//
