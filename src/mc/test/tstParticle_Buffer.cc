//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstParticle_Buffer.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 18:00:52 2001
 * \brief  Particle_Buffer testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "mc_test.hh"
#include "../Release.hh"
#include "../Particle_Buffer.hh"
#include "../Particle_Stack.hh"
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

using rtt_mc::Send_Particle_Buffer;
using rtt_mc::Recv_Particle_Buffer;
using rtt_mc::Particle_Buffer;
using rtt_mc::Particle_Containers;
using rtt_rng::Sprng;
using rtt_rng::Rnd_Control;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc_test::Dummy_Particle   DP;
typedef Particle_Containers<DP>::Bank Bank;

// random number seed
int seed = 392745;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void check_Dummy_Particle()
{
    vector<double> ref_rn(100);
    vector<char>   packed;
    {
	// make a Rnd_Control and get sprngs
	Rnd_Control control(seed);
	
	Sprng r1 = control.get_rn(10);
	Sprng rr = control.get_rn(10);
	
	// get reference random numbers
	for (int i = 0; i < ref_rn.size(); i++)
	    ref_rn[i] = rr.ran();

	// make a Dummy_Particle
	rtt_mc_test::Dummy_Particle p(101, 0.5, r1);

	// transport it
	p.transport(5, .2, 50);

	// pack it
	packed = p.pack();
    }

    // unpack a particle
    rtt_mc_test::Dummy_Particle up(packed);

    // check it
    if (!soft_equiv(up.get_wt(), 0.7)) ITFAILS;
    if (up.get_cell() != 106)          ITFAILS;

    Sprng r = up.get_rn();

    for (int i = 0; i < 20; i++)
	if (!soft_equiv(r.ran(), ref_rn[i+50])) ITFAILS;
    
    up.transport(0, 0.0, 10);

    // remember sprng is reference counted
    for (int i = 0; i < 20; i++)
	if (!soft_equiv(r.ran(), ref_rn[i+80])) ITFAILS;
	
    if (rtt_mc_test::passed)
	PASSMSG("Dummy_Particle ok.");
}

//---------------------------------------------------------------------------//

void setup_buffer()
{
    Rnd_Control control(seed);
    
    int before = control.get_num();
    int size_particle = rtt_mc_test::get_particle_size(control);
    if (control.get_num() != before) ITFAILS;

    // size the particle buffer
    Particle_Buffer<DP>::set_size_packed_particle(size_particle);

    if (rtt_mc_test::passed)
    {
	ostringstream message;
	message << "Particle_Buffer set with packed Dummy_Particle size of "
		<< Particle_Buffer<DP>::get_size_packed_particle();
	PASSMSG(message.str());
    }

    // try to make the buffer too big
    bool c = false;
    try
    {
	Particle_Buffer<DP>::set_maximum_num_particles(1000000000);
    }
    catch(rtt_dsxx::assertion &a)
    {
	ostringstream message;
	message << "Good, caught assertion, " << a.what();
	PASSMSG(message.str());
	c = true;
    }
    if (!c) 
    {
	FAILMSG("Failed to catch an illegal size assertion.");
    }
    if (Particle_Buffer<DP>::get_maximum_num_particles() != 1000) ITFAILS;
}

//---------------------------------------------------------------------------//

void test_buffer()
{
    // make a particle buffer
    Particle_Buffer<DP> buffer;
    if (!buffer.is_empty()) ITFAILS;

    vector<double> rnref1(100);
    vector<double> rnref2(100);

    Particle_Buffer<DP>::set_maximum_num_particles(3);

    {
	// make a rnd control
	Rnd_Control control(seed);

	Sprng r1  = control.get_rn(1);
	Sprng rr1 = control.get_rn(1);
	Sprng r2  = control.get_rn(2);
	Sprng rr2 = control.get_rn(2);
	Sprng r3(r1);

	for (int i = 0; i < 100; i++)
	{
	    rnref1[i] = rr1.ran();
	    rnref2[i] = rr2.ran();
	}
	
	// make some particles
	DP p1(1, 1.0, r1);
	DP p2(2, 2.0, r2);
	DP p3(3, 3.0, r3);

	// transport them
	p1.transport(0, 1.0, 25);
	p2.transport(0, 1.0, 25);
	p3.transport(0, 1.0, 10);

	// add them to the buffer
	buffer.buffer_particle(p1);
	buffer.buffer_particle(p2);
	buffer.buffer_particle(p3);
	
	if (buffer.get_num_particles_in_buffer() != 3) ITFAILS;
	if (!buffer.is_full())                         ITFAILS;

	// catch an illegal adding
	bool c = false;
	try
	{
	    buffer.buffer_particle(p1);
	}
	catch (const rtt_dsxx::assertion &caught)
	{
	    ostringstream message;
	    message << "Good, caught the following assertion, "
		    << caught.what();
	    PASSMSG(message.str());
	    c = true;
	}
	if (!c) 
	{
	    FAILMSG("Failed to catch an illegal buffer assertion.");
	}
    }

    // make a particle bank
    Particle_Containers<DP>::Bank bank;
    {
	// get particles out of the buffer and put them into the bank
	buffer.add_to_bank(bank);

	if (bank.size() != 3) ITFAILS;

	// check the particles
	SP<DP> p3 = bank.top();
	bank.pop();
	SP<DP> p2 = bank.top();
	bank.pop();
	SP<DP> p1 = bank.top();
	bank.pop();
	
	if (p1->get_cell() != 1) ITFAILS;
	if (p2->get_cell() != 2) ITFAILS;
	if (p3->get_cell() != 3) ITFAILS;

	if (p1->get_wt() != 2.0) ITFAILS;
	if (p2->get_wt() != 3.0) ITFAILS;
	if (p3->get_wt() != 4.0) ITFAILS;

	// transport the particles and check rn; remember even though 3
	// started out with a copy of r1, it got packed and unpacked, so at
	// this point it is its own number
	p1->transport(0, 0.0, 5);
	p2->transport(0, 0.0, 5);
	p3->transport(0, 0.0, 5);

	for (int i = 0; i < 60; i++)
	{
	    double r = p1->get_rn().ran();
	    if (!soft_equiv(r, rnref1[i+40])) ITFAILS;
	}

	for (int i = 0; i < 70; i++)
	{
	    double r = p2->get_rn().ran();
	    if (!soft_equiv(r, rnref2[i+30])) ITFAILS;
	}

	for (int i = 0; i < 60; i++)
	{
	    double r = p3->get_rn().ran();
	    if (!soft_equiv(r, rnref1[i+40])) ITFAILS;
	}
	
    }
    if (bank.size() != 0) ITFAILS;
	

    if (rtt_mc_test::passed)
	PASSMSG("Particle_Buffer buffering tests ok.");
}

//---------------------------------------------------------------------------//

void test_async()
{
    if (C4::nodes() != 4)
	return;
    
    // make a rnd control on each node
    Rnd_Control control(seed);

    // make reference random numbers
    vector<double> ref1(100);
    vector<double> ref2(100);
    vector<double> ref3(100);
    vector<double> ref4(100);
    vector<double> ref5(100);
    vector<double> ref6(100);

    if (Particle_Buffer<DP>::get_maximum_num_particles() != 3) ITFAILS;

    {
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);
	Sprng r3 = control.get_rn(3);
	Sprng r4 = control.get_rn(4);
	Sprng r5 = control.get_rn(5);
	Sprng r6 = control.get_rn(6);
	for (int i = 0; i < 100; i++)
	{
	    ref1[i] = r1.ran();
	    ref2[i] = r2.ran();
	    ref3[i] = r3.ran();
	    ref4[i] = r4.ran();
	    ref5[i] = r5.ran();
	    ref6[i] = r6.ran();
	}
    }

    // make a bank on each processor
    Bank bank;

    // make send/recv buffers for each processor
    // 0 -> 3
    // 1 -> 2 3
    // 2 -> 1 
    // 3 -> 0 2 
    
    // make send/recv buffers
    vector<Send_Particle_Buffer<DP> > send;
    vector<Recv_Particle_Buffer<DP> > recv;

    // make node maps
    vector<int> send_map;
    vector<int> recv_map;

    if (C4::node() == 0)
    {
	send.resize(1);
	recv.resize(1);
	send_map.resize(1);
	recv_map.resize(1);

	send_map[0] = 3;
	recv_map[0] = 3;

	// check that the buffer is not in use
	if (send[0].comm_status()) ITFAILS;
	if (recv[0].comm_status()) ITFAILS;

	// post a receive
	recv[0].post_arecv(recv_map[0]);

	// make a particle
	Sprng r = control.get_rn(1);
	DP particle(1, 1.0, r);

	// transport it
	particle.transport(5, .25, 15);

	// buffer it and send it
	send[0].buffer_particle(particle);
	send[0].post_asend(send_map[0]);
	if (!send[0].comm_status()) ITFAILS;
	if (send[0].is_empty())     ITFAILS;

	// try to catch an illegal adding assertion to send buffer
	bool caught = false;
	try
	{
	    send[0].buffer_particle(particle);
	}
	catch (const rtt_dsxx::assertion &ass)
	{
	    ostringstream message;
	    message << "Good, caught the following assertion, "
		    << ass.what() << " on processor " << C4::node();
	    PASSMSG(message.str());
	    caught = true;
	}
	if (!caught)
	{
	    FAILMSG("Failed to catch an illegal buffering assertion.");
	}

	// finish the send
	send[0].async_wait();
	if (send[0].comm_status())                      ITFAILS;
	if (!send[0].is_empty())                        ITFAILS;
	if (send[0].get_num_particles_in_buffer() != 0) ITFAILS;

	// wait on recv
	recv[0].async_wait();
	recv[0].add_to_bank(bank);
	if (!recv[0].is_empty()) ITFAILS;

	// check the particle
	if (bank.size() != 1) ITFAILS;

	SP<DP> rp = bank.top();
	bank.pop();

	if (rp->get_cell() != 5)            ITFAILS;
	if (!soft_equiv(rp->get_wt(), 2.5)) ITFAILS;

	for (int i = 0; i < 88; i++)
	{
	    double ran = rp->get_rn().ran();
	    double ref = ref1[i+12];
	    if (!soft_equiv(ran, ref)) ITFAILS;
	}
    }

    if (C4::node() == 1)
    {
	send.resize(2);
	recv.resize(1);
	send_map.resize(2);
	recv_map.resize(1);

	send_map[0] = 2;
	send_map[1] = 3;
	recv_map[0] = 2;

	// post a receive
	for (int i = 0; i < recv.size(); i++)
	    recv[i].post_arecv(recv_map[i]);

	// make a 2 particles
	{
	    Sprng r1 = control.get_rn(4);
	    Sprng r2 = control.get_rn(5);
	    
	    DP    p1(10, 1.2, r1);
	    DP    p2(12, 0.5, r2);
	    
	    // send (async) 
	    send[0].buffer_particle(p1);
	    send[1].buffer_particle(p2);
	}

	for (int i = 0; i < send.size(); i++)
	{
	    send[i].post_asend(send_map[i]);
	    send[i].async_wait();
	    if (!send[i].is_empty()) ITFAILS;
	}

	// receive the particles
	for (int i = 0; i < recv.size(); i++)
	{
	    recv[i].async_wait();
	    recv[i].add_to_bank(bank);
	}
	
	if (bank.size() != 1) ITFAILS;

	// check it 
	SP<DP> rp = bank.top();
	bank.pop();
	if (!bank.empty()) ITFAILS;

	if (rp->get_cell() != 113)          ITFAILS;
	if (!soft_equiv(rp->get_wt(), 1.1)) ITFAILS;

	for (int i = 0; i < 42; i++)
	{
	    double ran = rp->get_rn().ran();
	    double ref = ref6[i+58];
	    if (!soft_equiv(ran, ref)) ITFAILS;
	}
    }

    if (C4::node() == 2)
    {
	send.resize(1);
	recv.resize(2);
	send_map.resize(1);
	recv_map.resize(2);

	send_map[0] = 1;
	recv_map[0] = 1;
	recv_map[1] = 3;

	// post a receive
	for (int i = 0; i < recv.size(); i++)
	    recv[i].post_arecv(recv_map[i]);

	// send a particle to processor 1
	{
	    Sprng r = control.get_rn(6);
	    DP particle(1, 1.0, r);
	    particle.transport(112, .1, 58);
	    send[0].buffer_particle(particle);
	}

	for (int i = 0; i < send.size(); i++)
	{
	    send[i].post_asend(send_map[i]);
	    bool done = false;
	    while (!done)
	    {
		done = send[i].async_check();
	    }
	    if (!send[i].is_empty())    ITFAILS;
	    if (send[i].comm_status())  ITFAILS;
	}

	for (int i = 0; i < recv.size(); i++)
	{
	    if (!recv[i].comm_status()) ITFAILS;
	    bool done = false;
	    while (!done)
	    {
		done = recv[i].async_check();
	    }
	    if (recv[i].comm_status())  ITFAILS;
	    recv[i].add_to_bank(bank);
	    if (!recv[i].is_empty())    ITFAILS;
	}

	// check the particle
	if (bank.size() != 2) ITFAILS;

	while (!bank.empty())
	{
	    SP<DP> rp = bank.top();
	    bank.pop();

	    if (rp->get_cell() == 10)
	    {
		if (!soft_equiv(rp->get_wt(), 1.2)) ITFAILS;

		for (int i = 0; i < 100; i++)
		{
		    double ran = rp->get_rn().ran();
		    double ref = ref4[i];
		    if (!soft_equiv(ran, ref)) ITFAILS;
		}
	    }
	    
	    else if (rp->get_cell() == 7)
	    {
		if (!soft_equiv(rp->get_wt(), 3.25)) ITFAILS;

		for (int i = 0; i < 90; i++)
		{
		    double ran = rp->get_rn().ran();
		    double ref = ref2[i+10];
		    if (!soft_equiv(ran, ref)) ITFAILS;
		}
	    }

	    else
	    {
		FAILMSG("Processor 2 receives failed.");
	    }
	}
    }

    if (C4::node() == 3)
    {
	send.resize(2);
	recv.resize(2);
	send_map.resize(2);
	recv_map.resize(2);

	send_map[0] = 0;
	send_map[1] = 2;
	recv_map[0] = 0;
	recv_map[1] = 1;

	// post receives
	for (int i = 0; i < recv.size(); i++)
	    recv[i].post_arecv(recv_map[i]);

	// make two particles
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);
	DP    p1(2, 2.0, r1);
	DP    p2(3, 3.0, r2);

	// transport them
	p1.transport(3, .5, 12);
	p2.transport(4, .25, 10);

	// send out the buffers
	send[0].buffer_particle(p1);
	send[1].buffer_particle(p2);
	
	// blocking sends
	for (int i = 0; i < send.size(); i++)
	{
	    send[i].send_buffer(send_map[i]);
	    if (!send[i].is_empty()) ITFAILS;
	}

	// wait for receives
	for (int i = 0; i < recv_map.size(); i++)
	{
	    recv[i].async_wait();
	    recv[i].add_to_bank(bank);
	}

	// check that the bank has two particles
	if (bank.size() != 2) ITFAILS;

	// the first particle is from processor 0, next from 1; thus the
	// first particle out is from processor 1
	SP<DP> rp2 = bank.top();
	bank.pop();
	SP<DP> rp1 = bank.top();
	bank.pop();

	// check the particles
	if (rp1->get_cell() != 6)             ITFAILS;
	if (!soft_equiv(rp1->get_wt(), 1.25)) ITFAILS;

	for (int i = 0; i < 85; i++)
	{
	    double ran = rp1->get_rn().ran();
	    double ref = ref1[i+15];
	    if (!soft_equiv(ran, ref)) ITFAILS;
	}

	if (rp2->get_cell() != 12)             ITFAILS;
	if (!soft_equiv(rp2->get_wt(), 0.5))   ITFAILS;

	for (int i = 0; i < 100; i++)
	{
	    double ran = rp2->get_rn().ran();
	    double ref = ref5[i];
	    if (!soft_equiv(ran, ref)) ITFAILS;
	}
    }

    for (int i = 0; i < send.size(); i++)
    {
	if (send[i].comm_status()) ITFAILS;
	if (!send[i].is_empty())   ITFAILS;
    }

    for (int i = 0; i < recv.size(); i++)
    {
	if (recv[i].comm_status()) ITFAILS;
	if (!recv[i].is_empty())   ITFAILS;
    }

    if (rtt_mc_test::passed)
    {
	ostringstream message;
	message << "Async communication tests finished on processor "
		<< C4::node();
	PASSMSG(message.str());
    }
}

//---------------------------------------------------------------------------//

void test_async_again()
{
    if (C4::nodes() != 2) 
	return;
    
    // make a rnd control on each node
    Rnd_Control control(seed);

    // make reference random numbers
    vector<double> ref1(100);
    vector<double> ref2(100);

    if (Particle_Buffer<DP>::get_maximum_num_particles() != 3) ITFAILS;

    {
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);
	for (int i = 0; i < 100; i++)
	{
	    ref1[i] = r1.ran();
	    ref2[i] = r2.ran();
	}
    }

    // make a bank on each processor
    Bank bank;

    // make send/recv buffers
    Send_Particle_Buffer<DP> send;
    Recv_Particle_Buffer<DP> recv;
    
    int map = 1 - C4::node();

    // post receives
    recv.post_arecv(map);

    if (C4::node() == 0)
    {
	// make a particle
	{
	    Sprng r = control.get_rn(1);
	    DP    p(1, 1.5, r);
	    p.transport(2, .25, 10);
	    
	    // send it
	    send.buffer_particle(p);
	    send.post_asend(map);
	}

	// wait on the receives
	recv.async_wait();
	recv.add_to_bank(bank);

	// make sure the send is complete
	send.async_wait();

	// post a new receive
	recv.post_arecv(map);
	
	// wait on receive
	bool done = false;
	while (!done)
	    done = recv.async_check();
	recv.add_to_bank(bank);
	
	if (bank.size() != 2) ITFAILS;

	// check particles
	SP<DP> p1 = bank.top();
	bank.pop();
	SP<DP> p2 = bank.top();
	bank.pop();

	if (p1->get_cell() != 3)             ITFAILS;
	if (!soft_equiv(p1->get_wt(), 1.75)) ITFAILS;
	
	for (int i = 0; i < 90; i++)
	{
	    double r   = p1->get_rn().ran();
	    double ref = ref1[i+10];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}

	if (p2->get_cell() != 7)             ITFAILS;
	if (!soft_equiv(p2->get_wt(), 2.75)) ITFAILS;
	
	for (int i = 0; i < 90; i++)
	{
	    double r   = p2->get_rn().ran();
	    double ref = ref2[i+10];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }
    
    if (C4::node() == 1)
    {
	// make a Particle
	{
	    Sprng r = control.get_rn(2);
	    DP    p(4, 2.5, r);
	    p.transport(3, .25, 10);
	    
	    // send it
	    send.buffer_particle(p);
	    send.send_buffer(map);
	}

	// wait on receive
	recv.async_wait();
	recv.add_to_bank(bank);

	// pull the particle out of the bank and send it back
	send.buffer_particle(*bank.top());
	send.post_asend(map);
	send.async_wait();

	bank.pop();
	if (bank.size() != 0) ITFAILS;
    }

    if (recv.comm_status()) ITFAILS;
    if (send.comm_status()) ITFAILS;
    if (!recv.is_empty())   ITFAILS;
    if (!send.is_empty())   ITFAILS;

    if (rtt_mc_test::passed)
    {
	ostringstream message;
	message << "Multiple async communication tests finished on processor "
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
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// check out the dummy particle
	check_Dummy_Particle();

	// setup the buffer with the packed particle size
	setup_buffer();

	// test the buffer
	test_buffer();

	// test async communication
	test_async();
	test_async_again();
	
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstParticle_Buffer, " << ass.what()
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
	    cout << "**** tstParticle_Buffer Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstParticle_Buffer on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParticle_Buffer.cc
//---------------------------------------------------------------------------//
