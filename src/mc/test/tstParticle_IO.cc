//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstParticle_IO.cc
 * \author Thomas M. Evans
 * \date   Fri Jan  4 17:18:06 2002
 * \brief  Particle_IO test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../Particle_IO.hh"
#include "../Particle_Stack.hh"
#include "../Particle_Buffer.hh"
#include "../Release.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

using rtt_mc::Send_Particle_Buffer;
using rtt_mc::Recv_Particle_Buffer;
using rtt_mc::Particle_Buffer;
using rtt_mc::Particle_IO;
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

void file_output()
{
    // make a particle buffer
    Recv_Particle_Buffer<DP> buffer;
    if (!buffer.is_empty()) ITFAILS;

    // make a rnd control
    Rnd_Control control(seed);

    // set the size of the particle
    int size_particle = rtt_mc_test::get_particle_size(control);

    // size the particle buffer
    Particle_Buffer<DP>::set_size_packed_particle(size_particle);

    Particle_Buffer<DP>::set_maximum_num_particles(3);
    
    Sprng r1  = control.get_rn(1);
    Sprng rr1 = control.get_rn(1);
    Sprng r2  = control.get_rn(2);
    Sprng rr2 = control.get_rn(2);
    Sprng r3(r1);
	
    // make some particles
    DP p1(1, 1.0, r1);
    DP p2(2, 2.0, r2);
    DP p3(3, 3.0, r3);

    // transport them
    p1.transport(0, 1.0, 25);
    p2.transport(0, 1.0, 25);
    p3.transport(0, 1.0, 10);

    // add them to the buffer
    buffer.buffer_particle(p2);
    buffer.buffer_particle(p3);

    // make an output file
    ofstream file("particles");

    // write a particle to the file
    Particle_IO<DP>::write_particle(file, p1);
    
    // write the particle buffer to the file
    Particle_IO<DP>::write_Particle_Buffer(file, buffer);
}

//---------------------------------------------------------------------------//

void file_input()
{
    // reference numbers
    vector<double> rnref1(100);
    vector<double> rnref2(100);

    // make a rnd control
    Rnd_Control control(seed);

    {    
	Sprng rr1 = control.get_rn(1);
	Sprng rr2 = control.get_rn(2);

	for (int i = 0; i < 100; i++)
	{
	    rnref1[i] = rr1.ran();
	    rnref2[i] = rr2.ran();
	}
    }

    // make a bank
    Particle_Containers<DP>::Bank bank;

    // make the input file
    ifstream file("particles");

    // read particles from disk
    Particle_IO<DP>::read_particles(file, bank);

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

    // transport the particles and check rn; remember even though 3 started
    // out with a copy of r1, it got packed and unpacked, so at this point it
    // is its own number
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
    if (bank.size() != 0) ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Successfully recovered particles from input file.");
}

//---------------------------------------------------------------------------//

void check_empty()
{
    {
	ofstream ofile("empty_file");
	ifstream ifile("empty_file");

	// try to read the particles off the empty file

	// make a bank
	Particle_Containers<DP>::Bank bank;

	// read particles from disk
	Particle_IO<DP>::read_particles(ifile, bank);
 
	if (!bank.empty()) ITFAILS;
    }

    // check empty buffer write
    {
	ofstream ofile("empty_file");
	ifstream ifile("empty_file");

	Particle_Buffer<DP> buffer;

	Particle_IO<DP>::write_Particle_Buffer(ofile, buffer);

	// try to read the particles off the empty file

	// make a bank
	Particle_Containers<DP>::Bank bank;

	// read particles from disk
	Particle_IO<DP>::read_particles(ifile, bank);
 
	if (!bank.empty()) ITFAILS;
    }

    if (rtt_mc_test::passed)
	PASSMSG("Empty file read ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

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
	file_output();
	file_input();
	check_empty();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstParticle_IO, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstParticle_IO Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstParticle_IO on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParticle_IO.cc
//---------------------------------------------------------------------------//
