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
#include "../Topology.hh"
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
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::MC_Interface;
using dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

// Full replication topology
void replication_test()
{
    // get an OS_Mesh (2D 6 cells)
    SP<MC_Interface> interface(new MC_Interface());
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

    // make topology object
    Topology top(cpp, ppc, bc, "replication");

    // test num_cells
    if (top.num_cells() != mesh->num_cells()) ITFAILS;
    if (top.num_cells(C4::node()) != mesh->num_cells()) ITFAILS;
    
    // test num procs
    for (int cell = 1; cell <= top.num_cells(); cell++)
	if (top.num_procs(cell) != C4::nodes()) ITFAILS;

    // verify the parallel scheme
    if (top.get_parallel_scheme() != "replication") ITFAILS;

    // test global_cell functions
    for (int lc = 1; lc <= top.num_cells(C4::node()); lc++)
	if (top.global_cell(lc) != lc) ITFAILS;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int lc = 1; lc <= top.num_cells(proc); lc++)
	    if (top.global_cell(lc, proc) != lc) ITFAILS;

    // test local cell functions
    for (int gc = 1; gc <= top.num_cells(); gc++)
	if (top.local_cell(gc) != gc) ITFAILS;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int gc = 1; gc <= top.num_cells(); gc++)
	    if (top.local_cell(gc, proc) != gc) ITFAILS;

    // test get_cells function
    {
	vector<int> cell_list(mesh->num_cells());
	for (int i = 1; i <= cell_list.size(); i++)
	    cell_list[i-1] = i;

	// compare
	for (int np = 0; np < C4::nodes(); np++)
	{
	    vector<int> got_cells = top.get_cells(np);
	    if (got_cells != cell_list) ITFAILS;
	}
    }

    // test get_procs function
    {
	vector<int> proc_list(C4::nodes());
	for (int i = 0; i < C4::nodes(); i++)
	    proc_list[i] = i;

	// compare
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    vector<int> got_procs = top.get_procs(cell);
	    if (got_procs != proc_list) ITFAILS;
	}
    }
}

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

    // full replication test
    replication_test();

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
