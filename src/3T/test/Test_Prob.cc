//----------------------------------*-C++-*----------------------------------//
// Test_Prob.cc
// Geoffrey M. Furnish
// Wed Nov 19 16:18:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/test/Test_Prob.hh"
#include "3T/test/Test_3T.hh"
#include "3T/test/Run_DB.hh"
#include "3T/test/XYZ_Quadratic.hh"
#include "3T/test/XYZ_Trigonometric.hh"

#include "diffusion/Diffusion_DB.hh"

#include "mesh/Mesh_XYZ.hh"

#include "linalg/pcg_DB.hh"

#include "nml/Group.hh"

#include <string>
using std::string;

using dsxx::SP;

char *input_file = "test.in";
int verbose = 0;

void process_cli( int argc, char **argv )
{
    argc--, argv++;             // skip over program name.

    while( argc )
    {
        string cmd = *argv;

        cout << "processing option " << cmd << endl;

        if (cmd == "-f") {
            Assert( argc > 1 );
            input_file = argv[1];
            argc -= 2, argv += 2;
            continue;
        }

        if (cmd[0] == '-' && cmd.length() > 1 && cmd[1] == 'v') {
        // count the v's to determine the verbosity setting.
            for( int i=1; i < cmd.length(); i++ )
                if (cmd[i] == 'v') verbose++;
            argc--, argv++;
            continue;
        }

        if (verbose)
            cout << "unrecognized option: " << cmd << endl;

        argc--, argv++;
    }
}

SP<Test_Prob> Test_Prob_allocator( int argc, char *argv[] )
{
    process_cli( argc, argv );

    SP<Test_Prob> prob;

    NML_Group g( "test" );

    Run_DB rdb;
    rdb.setup_namelist( g );

    Diffusion_DB diffdb;
    diffdb.setup_namelist( g );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    Quad_Params qpdb;
    qpdb.setup_namelist( g );

    Trig_Params tpdb;
    tpdb.setup_namelist( g );

    pcg_DB pcg_db( "pcg" );
    pcg_db.setup_namelist( g );

    if (verbose)
        cout << "Reading input from " << input_file << ".\n";

    g.readgroup( input_file );
    g.writegroup( "test.out" );

    SP<Mesh_XYZ> spm = new Mesh_XYZ( mdb );

// Theoretically we could parse argc, argv to figure out which test problem
// to initiate.  For now, however, we just hardwire one.

    switch(rdb.test)
    {
    case Quad:
	prob = new Test_3T< Mesh_XYZ, XYZ_Quadratic >( spm, rdb, diffdb,
						       qpdb, pcg_db,
                                                       verbose );
	break;

    case Trig:
	prob = new Test_3T< Mesh_XYZ, XYZ_Trigonometric >( spm, rdb, diffdb,
							   tpdb, pcg_db,
                                                           verbose );
	break;

    default:
	throw "Unrecognized test problem specification.";
    }

    return prob;
}

//---------------------------------------------------------------------------//
//                              end of Test_Prob.cc
//---------------------------------------------------------------------------//
