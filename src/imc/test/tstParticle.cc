//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstParticle.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 14 17:07:23 2001
 * \brief  Particle and Particle_Buffer tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Particle.hh"
#include "../Particle_Buffer.hh"
#include "../Opacity.hh"
#include "../Tally.hh"
#include "../Release.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
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
using rtt_dsxx::SP;

// some global usefulness
int seed = 934223421;
Rnd_Control control(seed);
int proc;
int procs;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// easy test

template <typename POS>
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
    if (pt.status() != true)              ITFAILS;
    if (pt.get_ew() != 5.0)               ITFAILS;
    if (pt.get_descriptor() != POS::BORN) ITFAILS;
    if (pt.get_cell() != 100)             ITFAILS;
    if (pt.get_omega()[0] != 0.0)         ITFAILS;
    if (pt.get_omega()[1] != 1.0)         ITFAILS;
    if (pt.get_omega()[2] != 0.0)         ITFAILS;

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
    pt.set_descriptor(POS::SCATTER);
    pt.set_ew(11.0);
    pt.set_cell(30);
    POS pt2(r, omega, 11.0, 30, rnd2, 1, .11, POS::SCATTER);
    if (pt != pt2) ITFAILS;
}

//---------------------------------------------------------------------------//
// Pack test

template <typename POS>
void Particle_Pack()
{
    // typedef to nested Pack class so that static functions are usable
    typedef typename POS::Pack POS_Pack;

    // Make a particle
    vector<double> r(3, 0.0);
    vector<double> omega(3, 0.0);
    Sprng rnd = control.get_rn(10);
    r[1] = 2.5;
    omega[1] = 1.0;
    double ew = 5.0;
    int cell = 100;
    POS pt(r, omega, ew, cell, rnd);  // defaults for t_left and frac = 1.0

    // Check the status of Particle::Pack
    if (POS_Pack::get_int_size()!=0)    ITFAILS;
    if (POS_Pack::get_double_size()!=0) ITFAILS;
    if (POS_Pack::get_char_size()!=0)   ITFAILS;
    if (POS_Pack::get_setup())          ITFAILS;

    // Pack the particle

    // Get a smart pointer handle on the pack
    typename POS::SP_Pack packed(pt.pack());   

    // Check the setup and buffer sized of Particle::Pack
    if (POS_Pack::get_int_size()!=2)    ITFAILS;
    if (POS_Pack::get_double_size()!=9) ITFAILS;
    if (!POS_Pack::get_setup())         ITFAILS;

    // Check that the packed data is correct

    // double data
    const double *double_data=packed->double_begin();

    if (packed->get_double_size() != 9) ITFAILS;
    if (double_data[0] != 1.0)          ITFAILS;  // t_left default
    if (double_data[1] != ew )          ITFAILS;
    if (double_data[2] != 1.0)          ITFAILS;  // frac default
    if (double_data[3] != omega[0])     ITFAILS;  // omega[0]
    if (double_data[4] != omega[1])     ITFAILS;  // omega[1]
    if (double_data[5] != omega[2])     ITFAILS;  // omega[2]
    if (double_data[6] != r[0])         ITFAILS;  // r[0]
    if (double_data[7] != r[1])         ITFAILS;  // r[1]
    if (double_data[8] != r[2])         ITFAILS;  // r[2]

    // int data
    const int *int_data=packed->int_begin();

    if (packed->get_int_size() !=2)    ITFAILS;
    if (int_data[0] != cell)           ITFAILS; // cell
    if (int_data[1] != rnd.get_num())  ITFAILS; // random identfier

    // char data

    // Unpack the particle
    typename POS::SP_Particle unpacked(packed->unpack());

    // check data of unpacked particle
    if (unpacked->status() != true)              ITFAILS;
    if (unpacked->get_ew() != ew)                ITFAILS;
    if (unpacked->get_descriptor() != POS::BORN) ITFAILS;
    if (unpacked->get_cell() != cell)            ITFAILS;
    if (unpacked->get_omega()[0] != omega[0])    ITFAILS;
    if (unpacked->get_omega()[1] != omega[1])    ITFAILS;
    if (unpacked->get_omega()[2] != omega[2])    ITFAILS;
    if (unpacked->get_r()[0] != r[0])            ITFAILS;
    if (unpacked->get_r()[1] != r[1])            ITFAILS;
    if (unpacked->get_r()[2] != r[2])            ITFAILS;

    // Check to see that particles are identical
    if (*unpacked != pt) ITFAILS;
}

//---------------------------------------------------------------------------//
// Banking test

template <typename POS>
void Particle_Bank()
{

    // make a bank to put particles into
    typename Particle_Buffer<POS>::Bank bank;

    // make a bunch of particles
    for (int i = 1; i <= 5; i++)
    {
	vector<double> r(3, static_cast<double>(i));
	vector<double> omega(3, 0.0);
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
	vector<double> r(3, static_cast<double>(i));
	vector<double> omega(3, 0.0);
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

template <typename POS>
void Particle_Comm_Block()
{
    typedef Particle_Buffer<POS> PB;

    // make and test buffer
    int cs = control.get_size();
    //    PB buffer(8, 2, cs);
    PB buffer(3,control); // construct from dimension and random control object

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
	    vector<double> r(3, static_cast<double>(i));
	    vector<double> omega(3, 0.0);
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
	    if (send_comm[n-1].n_part != 5) ITFAILS;   // check for 5 particles
	    buffer.send_buffer(send_comm[n-1], n);
	    if (send_comm[n-1].n_part != 0) ITFAILS;   // check for empty buffer
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
	    vector<double> r(3, static_cast<double>(i));
	    vector<double> omega(3, 0.0);
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

template <typename POS>
void Particle_Comm_Async()
{
    typedef Particle_Buffer<POS> PB;

    // make particle buffers on each node
    //    PB buffer(8, 2, control.get_size());
    PB buffer(3, control);
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
	vector<double> r(3, static_cast<double>(proc));
	vector<double> omega(3);
	omega[2]         = 1.0;
	double ew        = (proc + 1) * 10;
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
    //    for (int i = 0; i < 10000 * i; i++); // doing nothing of consequence
    for (int i = 0; i < 10000  ; i++); // doing nothing of consequence

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
	vector<double> r(3, static_cast<double>(index[i]));
	vector<double> omega(3);
	omega[2]         = 1.0;
	double ew        = (index[i]+1) * 10;
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

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // proc definitions
    proc  = C4::node();
    procs = C4::nodes();

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

	// Particle tests 
	Particle_Basics<Particle<OS_Mesh> >();
	C4::gsync();
	
	Particle_Pack<Particle<OS_Mesh> >();
	C4::gsync();

	Particle_Bank<Particle<OS_Mesh> >();
	C4::gsync();
	
	Particle_Comm_Block<Particle<OS_Mesh> >();
	C4::gsync();
	
	Particle_Comm_Async<Particle<OS_Mesh> >();
	C4::gsync();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstParticle, " << ass.what()
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
	    cout << "**** tstParticle Test: PASSED on" 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstParticle on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParticle.cc
//---------------------------------------------------------------------------//
