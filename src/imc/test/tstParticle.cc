//----------------------------------*-C++-*----------------------------------//
// tstParticle.cc
// Thomas M. Evans
// Wed Apr 28 13:08:25 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Particle and Particle_Buffer Checks
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Particle.hh"
#include "../Particle_Buffer.hh"
#include "../Opacity.hh"
#include "../Tally.hh"
#include "../Release.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_imc::Particle;
using rtt_imc::Particle_Buffer;
using rtt_imc_test::IMC_Interface;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using dsxx::SP;

// some typedefs
typedef Particle<OS_Mesh> POS;
typedef Particle_Buffer<Particle<OS_Mesh> > PB;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_imc_test::fail(__LINE__);

// some global usefulness
int seed = 934223421;
Rnd_Control control(seed);
int proc;
int procs;

// Particle tests

//---------------------------------------------------------------------------//
// easy test

void Particle_Basics()
{
    // make a particle
    vector<double> r(3, 0.0);
    vector<double> omega(3, 0.0);
    Sprng rnd = control.get_rn(10);
    r[1] = 2.5;
    omega[1] = 1.0;
    double ew = 5.0;
    POS pt(r, omega, ew, 100, rnd);

    // tests
    if (pt.status() != true)   ITFAILS;
    if (pt.get_ew() != 5.0)    ITFAILS;
    if (pt.desc()   != "born") ITFAILS;
    if (pt.get_cell() != 100)  ITFAILS;

    // reset and check
    pt.kill_particle();
    if (pt.status() != false) ITFAILS;
    pt.reset_status();
    if (pt.status() != true)  ITFAILS;

    // make another particle to test equality
    Sprng rnd1 = control.get_rn(12);
    Sprng rnd2 = control.get_rn(12);
    pt.set_random(rnd1);
    pt.set_time_left(.11);
    pt.set_descriptor("test");
    pt.set_ew(11.0);
    pt.set_cell(30);
    POS pt2(r, omega, 11.0, 30, rnd2, 1, .11, "test");
    if (pt != pt2) ITFAILS;
}

//---------------------------------------------------------------------------//
// Banking test

void Particle_Bank()
{
    // make a bank to put particles into
    Particle_Buffer<POS>::Bank bank;

    // make a bunch of particles
    for (int i = 1; i <= 5; i++)
    {
	vector<double> r(3, i);
	vector<double> omega(3, 0);
	omega[0] = 1.0;
	double ew        = i * 10;
	double time_left = 1.0 / i;
	double fraction  = .5 / i;
	int cell         = i;
	Sprng rnd        = control.get_rn(i);
	SP<POS> particle(new POS(r, omega, ew, cell, rnd, fraction,
				 time_left)); 
	
	// place them into the bank
	bank.push(particle);
    }

    // check the bank
    if (bank.size() != 5) ITFAILS;

    for (int i = 5; i >= 1; i--)
    {
	vector<double> r(3, i);
	vector<double> omega(3, 0);
	omega[0] = 1.0;
	double ew        = i * 10;
	double time_left = 1.0 / i;
	double fraction  = .5 / i;
	int cell         = i;
	Sprng rnd        = control.get_rn(i);
	SP<POS> particle(new POS(r, omega, ew, cell, rnd, fraction,
				 time_left)); 
	if (*particle != *bank.top()) ITFAILS;
	bank.pop();
    }
    
    // check the bank, it should be empty
    if (bank.size())   ITFAILS;
    if (!bank.empty()) ITFAILS;
}

//---------------------------------------------------------------------------//
// pass some particles using a Particle_Buffer, BLOCKING COMMUNICATIONS

void Particle_Comm_Block()
{
    // make and test buffer
    int cs = control.get_size();
    PB buffer(8, 2, cs);

    if (PB::get_buffer_d() != 9000)      ITFAILS;
    if (PB::get_buffer_i() != 2000)      ITFAILS;
    if (PB::get_buffer_c() != cs * 1000) ITFAILS;

    PB::set_buffer_size(5);

    if (PB::get_buffer_d() != 45)     ITFAILS;
    if (PB::get_buffer_i() != 10)     ITFAILS;
    if (PB::get_buffer_c() != cs * 5) ITFAILS;

    // make and buffer five particles on processor 0
    if (!proc)
    {
	// we need to make a seperate Comm_Buffer for each node that we send
	// to because after a blocking send operation the data in the
	// Comm_Buffer is set to zero, so we need a Comm_Buffer for each node 
	// we send to, this will also test Comm_Buffer assignment
	PB::Comm_Vector send_comm(procs - 1);

	for (int i = 1; i <= 5; i++)
	{
	    vector<double> r(3, i);
	    vector<double> omega(3, 0);
	    omega[0] = 1.0;
	    double ew        = i * 10;
	    double time_left = 1.0 / i;
	    double fraction  = .5 / i;
	    int cell         = i;
	    Sprng rnd        = control.get_rn(i);
	    SP<POS> particle(new POS(r, omega, ew, cell, rnd, fraction,
				     time_left)); 
	
	    // place them into the send Comm_Buffer's
	    for (int j = 1; j < procs; j++)
		buffer.buffer_particle(send_comm[j-1], *particle);
	}

	// send out the comm buffer to ALL processors
	for (int n = 1; n < procs; n++)
	{
	    if (send_comm[n-1].n_part != 5) ITFAILS;
	    buffer.send_buffer(send_comm[n-1], n);
	    if (send_comm[n-1].n_part != 0) ITFAILS;
	}
    }

    // receive the buffers 
    if (proc)
    {
	PB::Bank recv_bank;
	SP<PB::Comm_Buffer> recv_comm = buffer.recv_buffer(0);
	if (recv_comm->n_part != 5) ITFAILS;

	buffer.add_to_bank(*recv_comm, recv_bank);
	if (recv_comm->n_part != 0) ITFAILS;

	// check the particles in the bank
	for (int i = 5; i >= 1; i--)
	{
	    vector<double> r(3, i);
	    vector<double> omega(3, 0);
	    omega[0] = 1.0;
	    double ew        = i * 10;
	    double time_left = 1.0 / i;
	    double fraction  = .5 / i;
	    int cell         = i;
	    Sprng rnd        = control.get_rn(i);
	    SP<POS> particle(new POS(r, omega, ew, cell, rnd, fraction,
				     time_left)); 
	
	    // check them
	    if (*particle != *recv_bank.top()) ITFAILS;
	    recv_bank.pop();
	}
	
	if (!recv_bank.empty()) ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// pass some particles using a Particle_Buffer, ASYNC COMMUNICATIONS

void Particle_Comm_Async()
{
    // make particle buffers on each node
    PB buffer(8, 2, control.get_size());
    PB::set_buffer_size(2);

    if (PB::get_buffer_d() != 18)                     ITFAILS;
    if (PB::get_buffer_i() != 4)                      ITFAILS;
    if (PB::get_buffer_c() != control.get_size() * 2) ITFAILS;
    
    // let's make a send/recv buffer on each node
    int send_buffers = procs - 1;
    int recv_buffers = procs - 1;
    PB::Comm_Vector recv_buf(recv_buffers);
    PB::Comm_Vector send_buf(send_buffers);

    // calculate the recv and send indices
    vector<int> index(procs - 1);
    {
	int indie = 0;
	for (int np = 0; np < procs; np++)
	    if (np != proc)
		index[indie++] = np;
    }
    Check (index.size() == recv_buf.size());
    Check (index.size() == send_buf.size());

    // post receives of the particles
    for (int i = 0; i < index.size(); i++)
    {
	if (recv_buf[i].n_part != 0) ITFAILS;
	buffer.post_arecv(recv_buf[i], index[i]);
    }

    // build two particles on each node
    {
	vector<double> r(3, proc);
	vector<double> omega(3);
	omega[2]         = 1.0;
	double ew        = proc * 10;
	double time_left = 1.0 / (proc + 1);
	double fraction  = .5 / (proc + 1);
	int cell         = proc + 1;
	Sprng r1         = control.get_rn(proc);
	Sprng r2         = control.get_rn(proc+1);
	int cell2        = proc + 2;
	POS p1(r, omega, ew, cell, r1, fraction, time_left);
	POS p2(r, omega, ew, cell2, r2, fraction, time_left);

	// pack the buffers
	for (int i = 0; i < index.size(); i++)
	{
	    buffer.buffer_particle(send_buf[i], p1);
	    buffer.buffer_particle(send_buf[i], p2);
	}
    }

    // make the processors do arbitrary amounts of work before sending its 
    // data
    for (int i = 0; i < 10000 * i; i++); // doing nothing of consequence

    // now lets send out the buffers
    for (int i = 0; i < index.size(); i++)
    {
	if (send_buf[i].n_part != 2) ITFAILS;
	buffer.asend_buffer(send_buf[i], index[i]);
	buffer.async_wait(send_buf[i]);
	if (send_buf[i].n_part != 0) ITFAILS;
    }

    // now lets check on our receive buffers until they are all in
    vector<bool> proc_in(index.size(), false);
    int all_in = 0;
    bool in = false;
    int iter;
    while (!in)
    {
	for (int i = 0; i < index.size(); i++)
	{
	    if (!proc_in[i])
		if (buffer.async_check(recv_buf[i]))
		{
		    proc_in[i] = true;
		    buffer.async_wait(recv_buf[i]);
		    if (recv_buf[i].n_part != 2) ITFAILS;
		    all_in++;
		}
	}
	if (all_in == index.size()) in = true;
	iter++;
    }

    // now we check the particles to make sure they are correct
    PB::Bank bank;
    for (int i = 0; i < index.size(); i++)
    {
	buffer.add_to_bank(recv_buf[i], bank);
	if (recv_buf[i].n_part != 0) ITFAILS;
    }
    if (bank.size() != 2 * index.size()) ITFAILS;

    // get particles from the bank and do our checks
    for (int i = index.size() - 1; i >= 0; i--)
    {
	// make reference particles
	vector<double> r(3, index[i]);
	vector<double> omega(3);
	omega[2]         = 1.0;
	double ew        = index[i] * 10;
	double time_left = 1.0 / (index[i] + 1);
	double fraction  = .5 / (index[i] + 1);
	int cell         = index[i] + 1;
	Sprng r1         = control.get_rn(index[i]);
	Sprng r2         = control.get_rn(index[i]+1);
	int cell2        = index[i] + 2;
	POS p1(r, omega, ew, cell, r1, fraction, time_left);
	POS p2(r, omega, ew, cell2, r2, fraction, time_left);

	if (*bank.top() != p2) ITFAILS;
	bank.pop();
	if (*bank.top() != p1) ITFAILS;
	bank.pop();
    }
    if (bank.size() != 0) ITFAILS;
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    // C4 Init
    C4::Init(argc, argv);
    proc  = C4::node();
    procs = C4::nodes();
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (proc == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl; 
	    C4::Finalize();
	    return 0;
	}

    // lets do the tests and catch assertions
    try
    {
	// Particle tests 
	Particle_Basics();
	C4::gsync();
	
	Particle_Bank();
	C4::gsync();
	
	Particle_Comm_Block();
	C4::gsync();
	
	Particle_Comm_Async();
	C4::gsync();
    }
    catch (const dsxx::assertion &ass)
    {
	cout << "Test: assertion failure at line " 
	     << ass.what() << endl;
	return 1;
    }
    catch(...)
    {
	cout << "HELP ME" << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "************************************" << endl;
    if (passed) 
    {
        cout << "**** Particle Self Test: PASSED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing Particle on node: " << proc << endl;

    // C4 end
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstParticle.cc
//---------------------------------------------------------------------------//
