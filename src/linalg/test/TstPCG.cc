//----------------------------------*-C++-*----------------------------------//
// TstPCG.cc
// Dave Nystrom
// Fri May  9 13:18:26 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include "../PCG_Ctrl.hh"
#include "../PCG_MatVec.hh"
#include "../PCG_PreCond.hh"

#include "test_utils.hh"
#include "TstPCG_MatVec.hh"
#include "TstPCG_PreCond.hh"

#include <iostream>
#include <string>

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    std::cout << progname << ": version " << version << std::endl;
}

//---------------------------------------------------------------------------//
// main
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
// Initialize C4.
    C4::Init( argc, argv );

// Provide for output of a version number
    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }

// Initialize some local variables
    int node  = C4::node();
    int nodes = C4::nodes();
    int ptype = C4::group();

// Now do the testing.
    int nxs = 12;
    int nys = 12;
    int nru = nxs * nys;

    using dsxx::SP;
    
    SP< PCG_MatVec<double> >  pcg_matvec(new TstPCG_MatVec<double>(nxs,nys));
    SP< PCG_PreCond<double> > pcg_precond(new TstPCG_PreCond<double>());

    PCG_Ctrl<double> pcg_ctrl(PCG_Ctrl<double>::GMRS);
    pcg_ctrl.setIparm(PCG_Ctrl<double>::ITSMAX, 1000);
    pcg_ctrl.setIparm(PCG_Ctrl<double>::NS2, 10);
    pcg_ctrl.setFparm(PCG_Ctrl<double>::ZETA, 0.001);
    pcg_ctrl.setFparm(PCG_Ctrl<double>::ALPHA, 0.1);
    pcg_ctrl.setOutputLevel(PCG_Ctrl<double>::LEVPRM);

    using dsxx::Mat1;
    
    Mat1<double> x(nru);
    Mat1<double> b(nru);

    double h = 1.0/(nxs+1);
    b = h*h;

    pcg_ctrl.solve( x, b, pcg_matvec, pcg_precond );

    // evaluate the results to see if it converged.
    
    Mat1<double> res(nru);

    // Get the results.
    
    pcg_matvec->MatVec(res, x);

    // The residual is the results minus the rhs.

    {
	C4::HTSyncSpinLock ht;

	const int ndigits = 2;
	
	for (int i = 0; i<nru; i++)
	{
	    using rtt_linalg::compare_reals;
	    using rtt_linalg::testMsg;
	    std::cout << res[i] << " " << b[i]
		      << " ---> "
		      << testMsg(compare_reals(res[i], b[i], ndigits))
		      << std::endl;
	}
    }

// Wrap up C4.
    C4::Finalize();

// Done.
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of TstPCG.cc
//---------------------------------------------------------------------------//
