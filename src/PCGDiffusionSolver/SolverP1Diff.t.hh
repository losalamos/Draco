//----------------------------------*-C++-*----------------------------------//
// SolverP1Diff.cc
// Randy M. Roberts
// Tue Sep 29 16:11:17 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "SolverP1Diff.hh"
#include "MatrixP1Diff.hh"
#include "MatVecP1Diff.hh"
#include "PreCondP1Diff.hh"
#include "ds++/Assert.hh"
#include <sstream>

namespace rtt_PCGDiffusionSolver
{

 template<class MT>
 SolverP1Diff<MT>::SolverP1Diff(const SP<const MT>& spMesh_,
				const pcg_DB& pcg_db)
     : spMesh(spMesh_), pcg_ctrl(pcgMethod(pcg_db.itmeth)),
     iqside(pcg_db.iqside)
 {
     pcg_ctrl.setOutputLevel(outputLevel(pcg_db.levout));
     pcg_ctrl.setStopTest(stopTest(pcg_db.ntest));

     // We should always want to use the incoming x as the initial guess
     Insist(pcg_db.iuinit == 1,
	    "Always use incoming x as the initial guess (pcg_db.iunit = 1).");

     pcg_ctrl.setUinit(uninit(pcg_db.iuinit));
     pcg_ctrl.setPrecon(precon(pcg_db.iqside));
     
     pcg_ctrl.setIparm(PCG_Ctrl::NOUT, pcg_db.nout);
     pcg_ctrl.setIparm(PCG_Ctrl::ITSMAX, pcg_db.itsmax);
     // pcg_ctrl.setIparm(PCG_Ctrl::MALLOC, pcg_db.malloc);
     // pcg_ctrl.setIparm(PCG_Ctrl::NWI, pcg_db.nwi);
     // pcg_ctrl.setIparm(PCG_Ctrl::NWF, pcg_db.nwf);
     // pcg_ctrl.setIparm(PCG_Ctrl::NEEDRC, pcg_db.needrc);
     pcg_ctrl.setIparm(PCG_Ctrl::NS1, pcg_db.ns1);
     pcg_ctrl.setIparm(PCG_Ctrl::NS2, pcg_db.ns2);

     // if (pcg_db.iuexac)
     // {
     //    .. set up the exact solution into uExact, to be tested
     //    pcg_ctrl.setUexact(const rtt_dsxx::Mat1<T>& uExact);
     // }

     pcg_ctrl.setLogical(PCG_Ctrl::ICKSTG, logical(pcg_db.ickstg));
     pcg_ctrl.setLogical(PCG_Ctrl::IDOT, logical(pcg_db.idot));
     pcg_ctrl.setLogical(PCG_Ctrl::ISTATS, logical(pcg_db.istats));

     pcg_ctrl.setFparm(PCG_Ctrl::CTIMER, pcg_db.ctimer);
     pcg_ctrl.setFparm(PCG_Ctrl::RTIMER, pcg_db.rtimer);
     pcg_ctrl.setFparm(PCG_Ctrl::FLOPSR, pcg_db.flopsr);
     pcg_ctrl.setFparm(PCG_Ctrl::ZETA, pcg_db.zeta);
     pcg_ctrl.setFparm(PCG_Ctrl::ALPHA, pcg_db.alpha);
 }

 template<class MT>
 void SolverP1Diff<MT>::solve(ccsf &phi, const SP<const Matrix> &spMatrix,
			      const ccsf &brhs)
 {
     SP<MatVec> spMatVec(new MatVec(spMatrix));
     
     SP<PreCond> spPreCond(new PreCond(spMatrix, iqside));

     // Now solve the matrix equation A.x = rhs.

     rtt_dsxx::Mat1<double> phitmp(phi.size());
     std::copy(phi.begin(), phi.end(), phitmp.begin());
    
     rtt_dsxx::Mat1<double> brhstmp(brhs.size());
     std::copy(brhs.begin(), brhs.end(), brhstmp.begin());

     rtt_dsxx::Mat1<double> &constbrhstmp = brhstmp;
     
     // pcg_ctrl.pcg_fe(phitmp, brhstmp, spMatVec, spPreCond);
     pcg_ctrl.solve(phitmp, constbrhstmp, spMatVec, spPreCond);
     
     std::copy(phitmp.begin(), phitmp.end(), phi.begin());
 }

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::Method
SolverP1Diff<MT>::pcgMethod(int itmeth)
{
    switch (itmeth)
    {
    case 1:
	return PCG_Ctrl::BAS;
    case 6:
	return PCG_Ctrl::GMRS;
    case 11:
	return PCG_Ctrl::CG;
    }

    std::ostringstream ost;
    ost << "SolverP1Diff: itmeth must be 1 -> BASIC, 6 -> GMRES, 11 -> CG"
	<< " found: " << itmeth;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::BAS;
}

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::OutputLevel
SolverP1Diff<MT>::outputLevel(int levout)
{
    switch (levout)
    {
    case PCG_Ctrl::LEV0:
	return PCG_Ctrl::LEV0;
    case PCG_Ctrl::LEVERR:
	return PCG_Ctrl::LEVERR;
    case PCG_Ctrl::LEVWRN:
	return PCG_Ctrl::LEVWRN;
    case PCG_Ctrl::LEVIT:
	return PCG_Ctrl::LEVIT;
    case PCG_Ctrl::LEVPRM:
	return PCG_Ctrl::LEVPRM;
    case PCG_Ctrl::LEVALG:
	return PCG_Ctrl::LEVALG;
    }
    std::ostringstream ost;
    ost << "SolverP1Diff: Invalid levout: " << levout;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::LEVALG;
}

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::StopTest
SolverP1Diff<MT>::stopTest(int ntest)
{
    switch (ntest)
    {
    case PCG_Ctrl::TSTUSR:
	return PCG_Ctrl::TSTUSR;
    case PCG_Ctrl::TSTEX:
	return PCG_Ctrl::TSTEX;
    case PCG_Ctrl::TSTDFA:
	return PCG_Ctrl::TSTDFA;
    case PCG_Ctrl::TST0:
	return PCG_Ctrl::TST0;
    case PCG_Ctrl::TSTSE:
	return PCG_Ctrl::TSTSE;
    case PCG_Ctrl::TSTSR:
	return PCG_Ctrl::TSTSR;
    case PCG_Ctrl::TSTSLR:
	return PCG_Ctrl::TSTSLR;
    case PCG_Ctrl::TSTSRR:
	return PCG_Ctrl::TSTSRR;
    case PCG_Ctrl::TSTRE:
	return PCG_Ctrl::TSTRE;
    case PCG_Ctrl::TSTRR:
	return PCG_Ctrl::TSTRR;
    case PCG_Ctrl::TSTRLR:
	return PCG_Ctrl::TSTRLR;
    case PCG_Ctrl::TSTRRR:
	return PCG_Ctrl::TSTRRR;
    }
    std::ostringstream ost;
    ost << "SolverP1Diff: Invalid ntest: " << ntest;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::TST0;
}

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::Uinit SolverP1Diff<MT>::uninit(int iuinit)
{
    switch (iuinit)
    {
    case PCG_Ctrl::USZERO:
	return PCG_Ctrl::USZERO;
    case PCG_Ctrl::UDFALT:
	return PCG_Ctrl::UDFALT;
    case PCG_Ctrl::UZERO:
	return PCG_Ctrl::UZERO;
    case PCG_Ctrl::UNZERO:
	return PCG_Ctrl::UNZERO;
    case PCG_Ctrl::USRAND:
	return PCG_Ctrl::USRAND;
    case PCG_Ctrl::UPRAND:
	return PCG_Ctrl::UPRAND;
    }
    std::ostringstream ost;
    ost << "SolverP1Diff: Invalid iuinit: " << iuinit;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::UNZERO;
}

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::Precon
SolverP1Diff<MT>::precon(int iqside)
{
    switch (iqside)
    {
    case PCG_Ctrl::QNONE:
	return PCG_Ctrl::QNONE;
    case PCG_Ctrl::QLEFT:
	return PCG_Ctrl::QLEFT;
    case PCG_Ctrl::QRIGHT:
	return PCG_Ctrl::QRIGHT;
    case PCG_Ctrl::QSPLIT:
	return PCG_Ctrl::QSPLIT;
    }
    std::ostringstream ost;
    ost << "SolverP1Diff: Invalid iqside: " << iqside;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::QNONE;
}

template<class MT>
typename SolverP1Diff<MT>::PCG_Ctrl::Logical
SolverP1Diff<MT>::logical(int value)
{
    switch (value)
    {
    case PCG_Ctrl::DFALT:
	return PCG_Ctrl::DFALT;
    case PCG_Ctrl::NO:
	return PCG_Ctrl::NO;
    case PCG_Ctrl::YES:
	return PCG_Ctrl::YES;
    }
    std::ostringstream ost;
    ost << "SolverP1Diff: Invalid logical: " << value;
    Insist(false, ost.str().c_str());
    // This return is just to shut up the compiler.
    return PCG_Ctrl::DFALT;
}

} // end namespace



//---------------------------------------------------------------------------//
//                              end of SolverP1Diff.cc
//---------------------------------------------------------------------------//
