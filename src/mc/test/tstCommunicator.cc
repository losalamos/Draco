//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstCommunicator.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:02:07 2001
 * \brief  Communicator and Comm_Builder tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "DD_Mesh.hh"
#include "../Particle_Buffer.hh"
#include "../Particle_Stack.hh"
#include "../Communicator.hh"
#include "../Communicator_Builder.hh"
#include "../Release.hh"
#include "../General_Topology.hh"
#include "../Rep_Topology.hh"
#include "../OS_Mesh.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_mc::Particle_Buffer;
using rtt_mc::Particle_Stack;
using rtt_mc::Communicator_Builder;
using rtt_mc::Communicator;
using rtt_mc::Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::General_Topology;
using rtt_mc::OS_Mesh;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc_test::Dummy_Particle PT;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// test Rep Communicator

void Rep_Comm_Check()
{
    // build a Rep Topology
    SP<Topology> topology(new Rep_Topology(9));

    // now build a communicator, the returned communicator should be null
    Communicator_Builder<PT> builder;
    SP<Communicator<PT> > comm = builder.build_Communicator(topology);
    if (comm) ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Replication communication test ok.");
}

//---------------------------------------------------------------------------//
// test DD Communicator and Build communicator

void DD_Comm_build()
{
    if (C4::nodes() != 4)
	return;

    // make a 4 processor/9 cell DD mesh
    SP<OS_Mesh> mesh = rtt_mc_test::build_Mesh();
    
    // make the topology that goes with it
    SP<Topology> topology = rtt_mc_test::build_Topology();

    // make a comm builder
    Communicator_Builder<PT> builder;

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

    if (rtt_mc_test::passed)
	PASSMSG("DD Communicator_Builder/Communicator check ok.");
}

// //---------------------------------------------------------------------------//
// // Transport some particles

void DD_Comm_transport()
{
    if (C4::nodes() != 4)
	return;

    // make a 4 processor/9 cell DD mesh
    SP<OS_Mesh> mesh = rtt_mc_test::build_Mesh();
    
    // make the topology that goes with it
    SP<Topology> topology = rtt_mc_test::build_Topology();

    // make a Random number controller
    SP<Rnd_Control> control(new Rnd_Control(349572));

    // reference numbers
    vector<double> ref10(100);
    vector<double> ref12(100);
    {
	Sprng r10 = control->get_rn(10);
	Sprng r12 = control->get_rn(12);

	for (int i = 0; i < 100; i++)
	{
	    ref10[i] = r10.ran();
	    ref12[i] = r12.ran();
	}
    }

    // set size of particle buffer
    int ps = rtt_mc_test::get_particle_size(*control);
    Particle_Buffer<PT>::set_size_packed_particle(ps);
    Particle_Buffer<PT>::set_maximum_num_particles(2);

    // make a comm builder
    Communicator_Builder<PT> builder;

    // make a communicator
    SP<Communicator<PT> > communicator =
	builder.build_Communicator(topology);

    if (communicator->get_send_size() != 0) ITFAILS;
    if (communicator->get_recv_size() != 0) ITFAILS;

    // post receives on all processors
    communicator->post();

    // make some particles on processor 0 headed for processor 1 and 2
    if (C4::node() == 0)
    {
	Sprng ran = control->get_rn(10);
	SP<PT> p1(new PT(1, 1.0, ran));
	SP<PT> p2(new PT(2, 2.0, ran));
	p2->transport(0,0.0,10);

	// make particles cross a boundary
	int indicator;
	p1->set_cell(-2);
	indicator = communicator->communicate(p1);
	if (indicator != -1)                    ITFAILS;
	if (communicator->get_send_size() != 1) ITFAILS;

	p2->set_cell(-1);
	indicator = communicator->communicate(p2);
	if (indicator != 1)                     ITFAILS;
	if (communicator->get_send_size() != 0) ITFAILS;
    }
    
    if (C4::node() == 1)
    {
	Particle_Stack<PT>::Bank bank;

	int arrived = 0;
	while (arrived != 3)
	    arrived += communicator->arecv_post(bank);
	if (arrived != 3) ITFAILS;

	if (bank.size() != 3)                   ITFAILS;
	if (communicator->get_recv_size() != 0) ITFAILS;

	// check particles
	while (!bank.empty())
	{
	    SP<PT> p = bank.top();
	    bank.pop();

	    // particle from proc 3
	    if (p->get_cell() == 2 && p->get_rn().get_num() == 12)
	    {
		if (!soft_equiv(p->get_wt(), 3.0)) ITFAILS;

		for (int i = 0; i < 80; i++)
		{
		    double r   = p->get_rn().ran();
		    double ref = ref12[i+20];
		    
		    if (!soft_equiv(r, ref)) ITFAILS;
		}
	    }

	    // first particle from proc 0
	    else if (p->get_cell() == 2 && p->get_rn().get_num() == 10)
	    {
		if (!soft_equiv(p->get_wt(), 1.0)) ITFAILS;

		for (int i = 0; i < 90; i++)
		{
		    double r   = p->get_rn().ran();
		    double ref = ref10[i+10];
		    
		    if (!soft_equiv(r, ref)) ITFAILS;
		}
	    }

	    // 2nd particle from proc 0
	    else if (p->get_cell() == 1)
	    {
		if (!soft_equiv(p->get_wt(), 2.0)) ITFAILS;

		for (int i = 0; i < 90; i++)
		{
		    double r   = p->get_rn().ran();
		    double ref = ref10[i+10];
		    
		    if (!soft_equiv(r, ref)) ITFAILS;
		}
	    }

	    else
	    {
		FAILMSG("Failed to get correct particle out of bank on 2");
	    }
	}
    }

    if (C4::node() == 3)
    {
	// this is potentially dangerous; however, because we don't receive
	// any particles on this node we can do this safely; also we are
	// setting the send size smaller than the value on other processors
	Particle_Buffer<PT>::set_maximum_num_particles(1);

	Sprng ran = control->get_rn(12);
	SP<PT> p1(new PT(1, 3.0, ran));
	p1->transport(0,0,20);
	
	int indicator;
	p1->set_cell(-1);
	indicator = communicator->communicate(p1);
	if (indicator != 1)                     ITFAILS;
	if (communicator->get_send_size() != 0) ITFAILS;
    }

    // sync everything up
    C4::gsync();

    if (!communicator->arecv_status()) ITFAILS;
    if (communicator->asend_status())  ITFAILS;
    if (communicator->get_send_size()) ITFAILS;
    if (communicator->get_recv_size()) ITFAILS;

    // end communication
    communicator->asend_end();
    communicator->arecv_end();

    if (communicator->arecv_status()) ITFAILS;
    if (communicator->asend_status()) ITFAILS;

    C4::gsync();

    if (rtt_mc_test::passed)
    {
	ostringstream message;
	message << "Communicator successfully communicated particles on "
		<< C4::node();
	PASSMSG(message.str());
    }
}

//---------------------------------------------------------------------------//
// check communicator functions on 4 proc

void DD_Comm_simple()
{
    if (C4::nodes() != 4)
	return;

    // make a 4 processor/9 cell DD mesh
    SP<OS_Mesh> mesh = rtt_mc_test::build_Mesh();
    
    // make the topology that goes with it
    SP<Topology> topology = rtt_mc_test::build_Topology();

    // make a Random number controller
    SP<Rnd_Control> control(new Rnd_Control(349572));

    // set size of particle buffer
    int ps = rtt_mc_test::get_particle_size(*control);
    Particle_Buffer<PT>::set_size_packed_particle(ps);
    Particle_Buffer<PT>::set_maximum_num_particles(2);

    // make a comm builder
    Communicator_Builder<PT> builder;

    // make a communicator
    SP<Communicator<PT> > communicator = builder.build_Communicator(topology);

    if (communicator->get_send_size() != 0) ITFAILS;
    if (communicator->get_recv_size() != 0) ITFAILS;

    // make a bank
    Particle_Stack<PT>::Bank bank;

    // >>> CHECK FUNCTIONS
    
    // post receives
    communicator->post();

    // make a particle on each node
    Sprng  ran = control->get_rn(1);
    SP<PT> particle(new PT(1, 1.0, ran)); 

    int received = 0;
    
    // add a particle to each buffer and then flush it
    if (C4::node() == 0)
    {
	// send to proc 2/local cell 1
	particle->set_cell(-3);

	if (communicator->communicate(particle) != -1) ITFAILS;
	
	vector<int> toproc = communicator->flush();
	if (toproc.size() != 1)            ITFAILS;
	if (toproc[0] != 2)                ITFAILS;
	if (communicator->get_send_size()) ITFAILS;
	if (communicator->asend_status())  ITFAILS;

	while (!received)
	    received = communicator->arecv_post(bank);

	if (bank.size() != 1) ITFAILS;
	if (!communicator->arecv_status()) ITFAILS;
	if (bank.top()->get_cell() != 2)   ITFAILS;

	if (rtt_mc_test::passed)
	    PASSMSG("Finished send/recv on proc 0.");
    }

    if (C4::node() == 1)
    {
	// send to proc 0/local cell 2
	particle->set_cell(-2);

	if (communicator->communicate(particle) != -1) ITFAILS;
	
	vector<int> toproc = communicator->flush();
	if (toproc.size() != 1)            ITFAILS;
	if (toproc[0] != 0)                ITFAILS;
	if (communicator->get_send_size()) ITFAILS;
	if (communicator->asend_status())  ITFAILS;

	while (!received)
	    received = communicator->arecv_post(bank);

	if (bank.size() != 1) ITFAILS;
	if (!communicator->arecv_status()) ITFAILS;
	if (bank.top()->get_cell() != 2)   ITFAILS;

	if (rtt_mc_test::passed)
	    PASSMSG("Finished send/recv on proc 1.");
    }
    
    if (C4::node() == 2)
    {
	// send to proc 3/local cell 2
	particle->set_cell(-4);

	if (communicator->communicate(particle) != -1) ITFAILS;
	
	vector<int> toproc = communicator->flush();
	if (toproc.size() != 1)            ITFAILS;
	if (toproc[0] != 3)                ITFAILS;
	if (communicator->get_send_size()) ITFAILS;
	if (communicator->asend_status())  ITFAILS;

	while (!received)
	    received = communicator->arecv_post(bank);

	if (bank.size() != 1) ITFAILS;
	if (!communicator->arecv_status()) ITFAILS;
	if (bank.top()->get_cell() != 1)   ITFAILS;

	if (rtt_mc_test::passed)
	    PASSMSG("Finished send/recv on proc 2.");
    }
    
    if (C4::node() == 3)
    {
	// send to proc 1/local cell 2
	particle->set_cell(-1);

	if (communicator->communicate(particle) != -1) ITFAILS;
	
	vector<int> toproc = communicator->flush();
	if (toproc.size() != 1)            ITFAILS;
	if (toproc[0] != 1)                ITFAILS;
	if (communicator->get_send_size()) ITFAILS;
	if (communicator->asend_status())  ITFAILS;

	while (!received)
	    received = communicator->arecv_post(bank);

	if (bank.size() != 1) ITFAILS;
	if (!communicator->arecv_status()) ITFAILS;
	if (bank.top()->get_cell() != 2)   ITFAILS;

	if (rtt_mc_test::passed)
	    PASSMSG("Finished send/recv on proc 3.");
    }

    C4::gsync();

    // end the communication (flush zero sized buffers and receive no new
    // particles) 
    communicator->flush_all();
    if (communicator->arecv_wait(bank)) ITFAILS;
    if (bank.size() != 1)               ITFAILS;

    if (communicator->get_send_size() != 0) ITFAILS;
    if (communicator->get_recv_size() != 0) ITFAILS;
    if (communicator->asend_status())       ITFAILS;
    if (communicator->arecv_status())       ITFAILS;

    C4::gsync();

    if (rtt_mc_test::passed)
    {
	ostringstream message;
	message << "Communicator simple test successfull on "
		<< C4::node();
	PASSMSG(message.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_mc::release() 
		 << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// rep test
	Rep_Comm_Check();

	// full DD test
	DD_Comm_build();
	DD_Comm_simple();
 	DD_Comm_transport();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstCommunicator, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstCommunicator Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstCommunicator on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstCommunicator.cc
//---------------------------------------------------------------------------//
