//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstCommunicator.cc
 * \author Thomas M. Evans
 * \date   Wed Jun  7 13:01:49 2000
 * \brief  Communicator and Comm_Builder test file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "DD_Mesh.hh"
#include "../Particle.hh"
#include "../Particle_Buffer.hh"
#include "../Communicator.hh"
#include "../Comm_Builder.hh"
#include "../Release.hh"
#include "mc/General_Topology.hh"
#include "mc/Rep_Topology.hh"
#include "mc/OS_Mesh.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_imc::Particle;
using rtt_imc::Comm_Builder;
using rtt_imc::Communicator;
using rtt_imc::Particle_Buffer;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::SP;

typedef Particle<OS_Mesh> PT;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// test Rep Communicator

void Rep_Comm_Check()
{
    // build a Rep Topology
    SP<Topology> topology(new Rep_Topology(9));

    // now build a communicator, the returned communicator should be null
    Comm_Builder<PT> builder;
    SP<Communicator<PT> > comm = builder.build_Communicator(topology);
    if (comm) ITFAILS;
}

//---------------------------------------------------------------------------//
// test DD Communicator and Build communicator

void DD_Comm_Check()
{
    if (C4::nodes() != 4)
	return;

    // make a 4 processor/9 cell DD mesh
    SP<OS_Mesh> mesh = rtt_imc_test::build_Mesh();
    
    // make the topology that goes with it
    SP<Topology> topology = rtt_imc_test::build_Topology();

    // make a comm builder
    Comm_Builder<PT> builder;

    // make a communicator
    SP<Communicator<PT> > communicator =
	builder.build_Communicator(topology);

    // some simple checks on the communicator
    if (C4::node() == 0)
    {
	if (communicator->num_send_nodes() != 2) ITFAILS;
	if (communicator->num_recv_nodes() != 2) ITFAILS;
    }
    else if (C4::node() == 1)
    {
	if (communicator->num_send_nodes() != 3) ITFAILS;
	if (communicator->num_recv_nodes() != 3) ITFAILS;
    }
    else if (C4::node() == 2)
    {
	if (communicator->num_send_nodes() != 3) ITFAILS;
	if (communicator->num_recv_nodes() != 3) ITFAILS;
    }
    if (C4::node() == 3)
    {
	if (communicator->num_send_nodes() != 2) ITFAILS;
	if (communicator->num_recv_nodes() != 2) ITFAILS;
    }

    // all communicators should have 0 particles in them
    if (communicator->get_send_size() != 0) ITFAILS;
    if (communicator->get_recv_size() != 0) ITFAILS;

    // check boundary cell size
    int bnd_cells = topology->get_boundary_cells(C4::node());
    if (communicator->num_bound_cells() != bnd_cells) ITFAILS;
    
    // define the boundary fields on each processor and test them
    vector<vector<int> > bnode(bnd_cells, vector<int>(1));
    vector<vector<int> > bcell(bnd_cells, vector<int>(1));
    vector<int>          com_nodes;
    if (C4::node() == 0)
    {
	bnode[0][0] = 1;
	bnode[1][0] = 1;
	bnode[2][0] = 2;
	
	bcell[0][0] = 1;
	bcell[1][0] = 2;
	bcell[2][0] = 1;
	
	if (bnode != communicator->get_b_node()) ITFAILS;
	if (bcell != communicator->get_b_cell()) ITFAILS;

	com_nodes.resize(2);
	com_nodes[0] = 1;
	com_nodes[1] = 2;
	
	// the send and recv nodes should be the same because this is full DD
	if (com_nodes != communicator->get_recv_nodes()) ITFAILS;
	if (com_nodes != communicator->get_send_nodes()) ITFAILS;
    }
    else if (C4::node() == 1)
    {
	bnode[0][0] = 0;
	bnode[1][0] = 0;
	bnode[2][0] = 2;
	bnode[3][0] = 2;
	bnode[4][0] = 3;
	
	bcell[0][0] = 1;
	bcell[1][0] = 2;
	bcell[2][0] = 1;
	bcell[3][0] = 2;
	bcell[4][0] = 1;

	if (bnode != communicator->get_b_node()) ITFAILS;
	if (bcell != communicator->get_b_cell()) ITFAILS;

	com_nodes.resize(3);
	com_nodes[0] = 0;
	com_nodes[1] = 2;
	com_nodes[2] = 3;
	
	// the send and recv nodes should be the same because this is full DD
	if (com_nodes != communicator->get_recv_nodes()) ITFAILS;
	if (com_nodes != communicator->get_send_nodes()) ITFAILS;
    }
    else if (C4::node() == 2)
    {
	bnode[0][0] = 0;
	bnode[1][0] = 1;
	bnode[2][0] = 1;
	bnode[3][0] = 3;
	bnode[4][0] = 3;
	
	bcell[0][0] = 2;
	bcell[1][0] = 1;
	bcell[2][0] = 2;
	bcell[3][0] = 2;
	bcell[4][0] = 3;

	if (bnode != communicator->get_b_node()) ITFAILS;
	if (bcell != communicator->get_b_cell()) ITFAILS;

	com_nodes.resize(3);
	com_nodes[0] = 0;
	com_nodes[1] = 1;
	com_nodes[2] = 3;
	
	// the send and recv nodes should be the same because this is full DD
	if (com_nodes != communicator->get_recv_nodes()) ITFAILS;
	if (com_nodes != communicator->get_send_nodes()) ITFAILS;
    }
    else if (C4::node() == 3)
    {
	bnode[0][0] = 1;
	bnode[1][0] = 2;
	bnode[2][0] = 2;
	
	bcell[0][0] = 2;
	bcell[1][0] = 1;
	bcell[2][0] = 2;

	if (bnode != communicator->get_b_node()) ITFAILS;
	if (bcell != communicator->get_b_cell()) ITFAILS;

	com_nodes.resize(2);
	com_nodes[0] = 1;
	com_nodes[1] = 2;
	
	// the send and recv nodes should be the same because this is full DD
	if (com_nodes != communicator->get_recv_nodes()) ITFAILS;
	if (com_nodes != communicator->get_send_nodes()) ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// Transport some particles

void DD_Comm()
{
    if (C4::nodes() != 4)
	return;

    // make a 4 processor/9 cell DD mesh
    SP<OS_Mesh> mesh = rtt_imc_test::build_Mesh();
    
    // make the topology that goes with it
    SP<Topology> topology = rtt_imc_test::build_Topology();

    // make a Random number controller
    SP<Rnd_Control> control(new Rnd_Control(349572));

    // make a Particle_Buffer
    SP<Particle_Buffer<PT> > buffer(new Particle_Buffer<PT>(*mesh,
							    *control)); 
    buffer->set_buffer_size(2);
    if (buffer->get_buffer_s() != 2) ITFAILS;

    // make a comm builder
    Comm_Builder<PT> builder;

    // make a communicator
    SP<Communicator<PT> > communicator =
	builder.build_Communicator(topology);

    if (communicator->get_send_size() != 0) ITFAILS;
    if (communicator->get_recv_size() != 0) ITFAILS;

    // post receives on all processors
    communicator->post(*buffer);

    // make some particles on processor 0 headed for processor 1 and 2
    if (C4::node() == 0)
    {
	vector<double> r(2, 0.0);
	vector<double> o(3, 0.0);
	Sprng ran = control->get_rn();
	SP<PT> p1(new PT(r, o, 1, 1, ran));
	SP<PT> p2(new PT(r, o, 2, 1, ran));

	// make particles cross a boundary
	int indicator;
	p1->set_cell(-2);
	indicator = communicator->communicate(*buffer, p1);
	if (indicator != -1)                    ITFAILS;
	if (communicator->get_send_size() != 1) ITFAILS;

	p2->set_cell(-1);
	indicator = communicator->communicate(*buffer, p2);
	if (indicator != 1)                     ITFAILS;
	if (communicator->get_send_size() != 0) ITFAILS;
    }
    
    if (C4::node() == 1)
    {
	Particle_Buffer<PT>::Bank bank;

	int arrived = 0;
	while (arrived != 2)
	    arrived += communicator->arecv_post(*buffer, bank);
	if (arrived != 2) ITFAILS;

	if (bank.size() != 3)                   ITFAILS;
	if (communicator->get_recv_size() != 0) ITFAILS;
    }

    if (C4::node() == 3)
    {
	buffer->set_buffer_size(1);

	vector<double> r(2, 0.0);
	vector<double> o(3, 0.0);
	Sprng ran = control->get_rn();
	SP<PT> p1(new PT(r, o, 1, 1, ran));
	
	int indicator;
	p1->set_cell(-1);
	indicator = communicator->communicate(*buffer, p1);
	if (indicator != 1)                     ITFAILS;
	if (communicator->get_send_size() != 0) ITFAILS;
    }

    // sync everything up
    C4::gsync();

    if (!communicator->arecv_status(*buffer)) ITFAILS;
    if (communicator->asend_status(*buffer))  ITFAILS;
    if (communicator->get_send_size())        ITFAILS;
    if (communicator->get_recv_size())        ITFAILS;

    // end communication
    communicator->asend_end(*buffer);
    communicator->arecv_end(*buffer);

    if (communicator->arecv_status(*buffer)) ITFAILS;
    if (communicator->asend_status(*buffer)) ITFAILS;
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl; 
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// rep test
	Rep_Comm_Check();

	// full DD test
	DD_Comm_Check();
	DD_Comm();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "Test: assertion failure at line " 
	     << ass.what() << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "****************************************" << endl; 
    if (passed) 
    {
        cout << "**** Communicator Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "****************************************" << endl;
    cout << endl;

    cout << "Done testing Communicator on node: " << C4::node() << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstCommunicator.cc
//---------------------------------------------------------------------------//
