//----------------------------------*-C++-*----------------------------------//
// TstPCG.cc
// Dave Nystrom
// Fri May  9 13:18:26 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/global.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

//#include "base/Log.hh"

#include "../pcg_DB.hh"
#include "../PCG_Ctrl.hh"
#include "../PCG_MatVec.hh"
#include "../PCG_PreCond.hh"

#include "tstpcg_DB.hh"
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

// Now, try to read in some namelist input.
    NML_Group *g;

    tstpcg_DB tstpcg_db;
    pcg_DB    pcg_db( "pcg" );

    g = new NML_Group( "TstPCG.in" );

    tstpcg_db.setup_namelist( *g );
    pcg_db   .setup_namelist( *g );

    g->readgroup ( "TstPCG.in"  );
    g->writegroup( "TstPCG.out" );

// Now, try writing to a log file.
//     Log tstlog;
//     tstlog.init( "TstPCG.log" );
//     tstlog << "*** Log File for TstPCG ***\n";

// Now do the testing.
    int nxs = tstpcg_db.nxs;
    int nys = tstpcg_db.nys;
    int nru = nxs * nys;

    using dsxx::SP;
    
    SP< PCG_MatVec<double> >  pcg_matvec  = new TstPCG_MatVec<double>(nxs,nys);
    SP< PCG_PreCond<double> > pcg_precond = new TstPCG_PreCond<double>();

    PCG_Ctrl<double> pcg_ctrl( pcg_db, nru );

    using dsxx::Mat1;
    
    Mat1<double> x(nru);
    Mat1<double> b(nru);

    double h = 1.0/(nxs+1);
    b = h*h;

    pcg_ctrl.pcg_fe( x, b, pcg_matvec, pcg_precond );

// Wrap up C4.
    C4::Finalize();

// Done.
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of TstPCG.cc
//---------------------------------------------------------------------------//
