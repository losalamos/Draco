//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstParticle.cc
 * \author Thomas M. Evans
 * \date   Mon Jan  7 18:41:29 2002
 * \brief  Particle test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "../Release.hh"
#include "../Particle.hh"
#include "mc/OS_Mesh.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc::Particle;
using rtt_mc::OS_Mesh;
using rtt_rng::Sprng;
using rtt_rng::Rnd_Control;
using rtt_dsxx::soft_equiv;

typedef Particle<OS_Mesh> PT;

// random number seed
int seed = 395731;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void simple_test()
{
    // random controller
    Rnd_Control control(seed);
    
    // make reference random numbers
    vector<double> ref1(100);
    vector<double> ref2(100);
    vector<double> ref3(100);
    {
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);
	Sprng r3 = control.get_rn(3);

	for (int i = 0; i < 100; i++)
	{
	    ref1[i] = r1.ran();
	    ref2[i] = r2.ran();
	    ref3[i] = r3.ran();
	}
    }

    // make a particle
    vector<double> r(3);
    vector<double> o(3);
    double         ew   = 0.45;
    double         tl   = 0.15;
    double         fr   = 0.5;
    int            cell = 10;
    Sprng          rnd  = control.get_rn(1);

    // throw some random numbers
    for (int i = 0; i < 15; i++)
	rnd.ran();

    r[0] = 0.5;
    r[1] = 0.75;
    r[2] = 4.0;
    o[0] = 0.5;
    o[1] = 0.5;
    o[2] = std::sqrt(2.0)/2.0;
    PT pref(r, o, ew, cell, rnd, fr, tl);

    // check the particle after a copy
    {
	PT p = pref;

	// check it
	vector<double> rf(3);
	vector<double> of(3);

	rf[0] = 0.5;
	rf[1] = 0.75;
	rf[2] = 4.0;
	of[0] = 0.5;
	of[1] = 0.5;
	of[2] = std::sqrt(2.0)/2.0;
	
	vector<double> r = p.get_r();
	vector<double> o = p.get_omega();
	
	if (!soft_equiv(r.begin(), r.end(), rf.begin(), rf.end())) ITFAILS;
	if (!soft_equiv(o.begin(), o.end(), of.begin(), of.end())) ITFAILS;
	if (p.get_cell() != 10)                                    ITFAILS;
	if (!soft_equiv(p.get_ew(), 0.45))                         ITFAILS;
	if (!p.status())                                           ITFAILS;
	if (p.get_descriptor() != PT::BORN)                        ITFAILS;

	// check random numbers; the stream is the same as the original
	// particle 
	for (int i = 0; i < 20; i++)
	{
	    double r   = p.get_random().ran();
	    double ref = ref1[i+15];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }

    // check that the original particle is on the same stream
    for (int i = 0; i < 65; i++)
    {
	double r   = pref.get_random().ran();
	double ref = ref1[i+35];
	if (!soft_equiv(r, ref)) ITFAILS;
    }

    // check set functions
    
    

    if (rtt_imc_test::passed)
	PASSMSG("Particle simple tests ok.");
}

//---------------------------------------------------------------------------//

void pack_test()
{
    // packed particle
    vector<char> packed;
    
    // random controller
    Rnd_Control control(seed);
    
    // make reference random numbers
    vector<double> ref1(100);
    {
	Sprng r1 = control.get_rn(1);

	for (int i = 0; i < 100; i++)
	    ref1[i] = r1.ran();
    }

    // pack a particle
    {
	vector<double> r(3);
	vector<double> o(3);
	double         ew   = 0.45;
	double         tl   = 0.15;
	double         fr   = 0.5;
	int            cell = 10;
	Sprng          rnd  = control.get_rn(1);

	// throw some random numbers
	for (int i = 0; i < 15; i++)
	    rnd.ran();

	r[0] = 0.5;
	r[1] = 0.75;
	r[2] = 4.0;
	o[0] = 0.5;
	o[1] = 0.5;
	o[2] = std::sqrt(2.0)/2.0;
	PT p(r, o, ew, cell, rnd, fr, tl);

	// pack the particle
	packed = p.pack();

	// throw more random numbers
	rnd.ran();
    }

    // unpack it and test
    {
	PT p(packed);

	// check it
	vector<double> rf(3);
	vector<double> of(3);

	rf[0] = 0.5;
	rf[1] = 0.75;
	rf[2] = 4.0;
	of[0] = 0.5;
	of[1] = 0.5;
	of[2] = std::sqrt(2.0)/2.0;
	
	vector<double> r = p.get_r();
	vector<double> o = p.get_omega();
	
	if (!soft_equiv(r.begin(), r.end(), rf.begin(), rf.end())) ITFAILS;
	if (!soft_equiv(o.begin(), o.end(), of.begin(), of.end())) ITFAILS;
	if (p.get_cell() != 10)                                    ITFAILS;
	if (!soft_equiv(p.get_ew(), 0.45))                         ITFAILS;
	if (!p.status())                                           ITFAILS;
	if (p.get_descriptor() != PT::UNPACKED)                    ITFAILS;

	// check random numbers; the stream continues at the point it was
	// packed 
	for (int i = 0; i < 85; i++)
	{
	    double r   = p.get_random().ran();
	    double ref = ref1[i+15];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }

    if (rtt_imc_test::passed)
	PASSMSG("Particle pack test ok.");
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
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	simple_test();
	pack_test();
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
	    cout << "**** tstParticle Test: PASSED on " 
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
