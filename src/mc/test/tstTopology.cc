//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstTopology.cc
 * \author Thomas M. Evans
 * \date   Thu Nov 18 16:02:11 1999
 * \brief  Topology test code.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../General_Topology.hh"
#include "../Rep_Topology.hh"
#include "../OS_Mesh.hh"
#include "../OS_Builder.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

using rtt_mc::Topology;
using rtt_mc::General_Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::Parser;
using rtt_dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__, __FILE__);

//---------------------------------------------------------------------------//
// REPLICATION TOPOLOGY TESTS
//---------------------------------------------------------------------------//
// full replication/General_Topology test

void test_Replication()
{
    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();

    // build a topology based on full replication
    Topology::vf_int cpp(C4::nodes());
    Topology::vf_int ppc(mesh->num_cells());
    Topology::vf_int bc;

    // fill up cpp vector
    for (int np = 0; np < cpp.size(); np++)
    {
	cpp[np].resize(mesh->num_cells());
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	    cpp[np][cell-1] = cell;
    }

    // fill up ppc vector
    for (int cell = 1; cell <= ppc.size(); cell++)
    {
	ppc[cell-1].resize(C4::nodes());
	for (int np = 0; np < C4::nodes(); np++)
	    ppc[cell-1][np] = np;
    }

    // build and test General_Topology for replication
    General_Topology general(cpp, ppc, bc, "replication");
    if ( !rtt_mc_test::topology_replication_test(mesh, general) ) ITFAILS;

    // build and test a SP<Topology> for replication
    SP<Topology> spgen(new General_Topology(cpp, ppc, bc, "replication"));
    if ( !rtt_mc_test::topology_replication_test(mesh, *spgen) ) ITFAILS;

    // build and test Rep_Topology
    Rep_Topology rep(mesh->num_cells());
    if ( !rtt_mc_test::topology_replication_test(mesh, rep) ) ITFAILS;

    // build and test SP<Topology> for replication
    SP<Topology> sprep(new Rep_Topology(mesh->num_cells()));
    if ( !rtt_mc_test::topology_replication_test(mesh, *sprep) ) ITFAILS;
}

//---------------------------------------------------------------------------//
// DD TESTS
//---------------------------------------------------------------------------//

void test_DD()
{
    // only perform this test on two processors
    if (C4::nodes() != 2) 
	return;

    // get an OS_Mesh (2D 6 cells)
    SP<Parser> interface(new Parser());
    OS_Builder builder(interface);
    SP<OS_Mesh> mesh = builder.build_Mesh();

    // build a topology based on full DD
    Topology::vf_int cpp(C4::nodes());
    Topology::vf_int ppc(mesh->num_cells());
    Topology::vf_int bc(C4::nodes());

    // fill up cells per processor array
    cpp[0].resize(3);
    cpp[1].resize(3);

    cpp[0][0] = 1;
    cpp[0][1] = 2;
    cpp[0][2] = 3;
    cpp[1][0] = 4;
    cpp[1][1] = 5;
    cpp[1][2] = 6;

    // fill up processor per cells array
    for (int i = 0; i < ppc.size(); i++)
	ppc[i].resize(1);

    ppc[0][0] = 0;
    ppc[1][0] = 0;
    ppc[2][0] = 0;
    ppc[3][0] = 1;
    ppc[4][0] = 1;
    ppc[5][0] = 1;

    // fill up boundary cells array
    bc[0].resize(3);
    bc[1].resize(3);
    
    bc[0][0] = 4;
    bc[0][1] = 5;
    bc[0][2] = 6;   
    bc[1][0] = 1;
    bc[1][1] = 2;
    bc[1][2] = 3;

    // build topology
    General_Topology topology(cpp, ppc, bc, "DD");

    // test topology
    if ( !rtt_mc_test::topology_DD_test(mesh, topology) ) ITFAILS;
}

//---------------------------------------------------------------------------//
// MAIN
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (!C4::node())
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    // full replication tests of General_Topology and Rep_Topology
    test_Replication();

    // DD tests of General_Topology
    test_DD();

    // status of test
    cout << endl;
    cout <<     "************************************" << endl;
    if (passed) 
    {
        cout << "**** Topology Self Test: PASSED on " << C4::node()
	     << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing Topology on node " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstTopology.cc
//---------------------------------------------------------------------------//
