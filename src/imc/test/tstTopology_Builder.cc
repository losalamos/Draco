//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstTopology_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:12:55 2001
 * \brief  Topology_Builder test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Topology_Builder.hh"
#include "../Release.hh"
#include "mc/Topology.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>

using namespace std;

using rtt_imc_test::IMC_Flat_Interface;
using rtt_imc_test::Parser;
using rtt_imc::Topology_Builder;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// build and test a full replication Topology built by Topology_Builder

void replication_test()
{	    
    // build a mesh
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // get the dummy interface with a capacity of 6 cells (the full mesh)
    SP<IMC_Flat_Interface> interface(new IMC_Flat_Interface(mb, 6));
    
    // build the Topology Builder and full replication topology
    Topology_Builder<OS_Mesh> tb(interface);
    SP<Topology> topology = tb.build_Topology(mesh);
    
    // check to make sure we have a Rep_Topology instance of Topology
    if (typeid(*topology) != typeid(Rep_Topology))   ITFAILS;
    if (typeid(*topology) == typeid(Topology))       ITFAILS;
    if (typeid(topology.bp()) != typeid(Topology *)) ITFAILS;
    
    // check to make sure the Topology is full replication
    if (topology->get_parallel_scheme() != "replication") ITFAILS;
    
    // check the number of cells
    if (topology->num_cells() != mesh->num_cells()) ITFAILS;
    
    // we have already checked the rest of Rep_Topology thoroughly in
    // mc/test/tstTopology so we will not repeat it here-->we can do this
    // because Rep_Topology has no data associated with it besides the number 
    // of global cells, this since we checked the functions previously we
    // know that this thing works
}

//---------------------------------------------------------------------------//
// build and test a DD Topology built by Topology builder

void DD_test()
{
    // this test is only for 2 processor tests
    if (C4::nodes() != 2) 
	return;

    // build a mesh
    SP<Parser> parser(new Parser("OS_Input"));
    SP<OS_Builder> mb(new OS_Builder(parser));
    SP<OS_Mesh> mesh = mb->build_Mesh();

    // get the dummy interface with a capacity of 3 cells (2 processor)
    SP<IMC_Flat_Interface> interface(new IMC_Flat_Interface(mb, 3));

    // build the Topology builder and full replication topology
    Topology_Builder<OS_Mesh> tb(interface);
    SP<Topology> topology = tb.build_Topology(mesh);

    // check to make sure we have a General_Topology instance of Topology
    if (typeid(*topology) != typeid(General_Topology)) ITFAILS;
    if (typeid(topology.bp()) != typeid(Topology *))   ITFAILS;

    // check to make sure Topology is full DD
    if (topology->get_parallel_scheme() != "DD") ITFAILS;

    // check to make sure the full DD topology is properly constructed, we
    // only check the topology on node 0 because that is what
    // Topology_Builder does->to check Topology functionality on other nodes
    // see the mc/test topology test

    Topology &top = *topology;
    
    // check local cells
    if (top.local_cell(1, 0) != 1) ITFAILS;
    if (top.local_cell(2, 0) != 2) ITFAILS;
    if (top.local_cell(3, 0) != 3) ITFAILS;
    if (top.local_cell(4, 0) != 0) ITFAILS;
    if (top.local_cell(5, 0) != 0) ITFAILS;
    if (top.local_cell(6, 0) != 0) ITFAILS;
    if (top.local_cell(1, 1) != 0) ITFAILS;
    if (top.local_cell(2, 1) != 0) ITFAILS;
    if (top.local_cell(3, 1) != 0) ITFAILS;
    if (top.local_cell(4, 1) != 1) ITFAILS;
    if (top.local_cell(5, 1) != 2) ITFAILS;
    if (top.local_cell(6, 1) != 3) ITFAILS;

    // check global cells
    if (top.num_cells(0) != 3) ITFAILS;
    if (top.num_cells(1) != 3) ITFAILS;
    
    if (top.global_cell(1, 0) != 1) ITFAILS;
    if (top.global_cell(2, 0) != 2) ITFAILS;
    if (top.global_cell(3, 0) != 3) ITFAILS;
    if (top.global_cell(1, 1) != 4) ITFAILS;
    if (top.global_cell(2, 1) != 5) ITFAILS;
    if (top.global_cell(3, 1) != 6) ITFAILS;

    // check cell lists
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (top.num_procs(cell) != 1) ITFAILS;

    vector<int> cell_list_0(3);
    vector<int> cell_list_1(3);
    vector<int> proc_list_0(1, 0);
    vector<int> proc_list_1(1, 1);

    cell_list_0[0] = 1;
    cell_list_0[1] = 2;
    cell_list_0[2] = 3;
    cell_list_1[0] = 4;
    cell_list_1[1] = 5;
    cell_list_1[2] = 6;

    if (top.get_cells(0) != cell_list_0) ITFAILS;
    if (top.get_cells(1) != cell_list_1) ITFAILS;
    
    if (top.get_procs(1) != proc_list_0) ITFAILS;
    if (top.get_procs(2) != proc_list_0) ITFAILS;
    if (top.get_procs(3) != proc_list_0) ITFAILS;
    if (top.get_procs(4) != proc_list_1) ITFAILS;
    if (top.get_procs(5) != proc_list_1) ITFAILS;
    if (top.get_procs(6) != proc_list_1) ITFAILS;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_imc::release() 
		 << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// Topology_Builder is a serial node object builder
	if (!C4::node())
	{
	    // full replication test
	    replication_test();
	
	    // full DD test
	    DD_test();
	}
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstTopology_Builder, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstTopology_Builder Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstTopology_Builder on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstTopology_Builder.cc
//---------------------------------------------------------------------------//
