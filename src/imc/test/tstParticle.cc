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

#include "IMC_Test.hh"
#include "imc_test.hh"
#include "../Release.hh"
#include "../Particle.hh"
#include "../Opacity.hh"
#include "../Tally.hh"
#include "../Global.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
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
using rtt_imc::Opacity;
using rtt_imc::Tally;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_rng::Sprng;
using rtt_rng::Rnd_Control;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::SP;

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
    {
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);

	for (int i = 0; i < 100; i++)
	{
	    ref1[i] = r1.ran();
	    ref2[i] = r2.ran();
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

	if (p != pref) ITFAILS;

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
    {
	Sprng r2 = control.get_rn(2);

	pref.set_random(r2);

	for (int i = 0; i < 50; i++)
	{
	    double r   = pref.get_random().ran();
	    double ref = ref2[i];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }

    {
	pref.set_cell(12);
	pref.set_ew(0.25);

	if (pref.get_cell() != 12)               ITFAILS;
	if (!soft_equiv(pref.get_ew(), 0.25))    ITFAILS;

	pref.kill_particle();

	if (pref.status())                       ITFAILS;

	pref.reset_status();

	if (!pref.status())                      ITFAILS;

	if (pref.get_descriptor() != PT::BORN)   ITFAILS;
	pref.set_descriptor(PT::CENSUS);
	if (pref.get_descriptor() != PT::CENSUS) ITFAILS;
    }
    
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

	// check the size packed state function
	if (PT::get_packed_particle_size(2, control) == packed.size()) ITFAILS;
	if (PT::get_packed_particle_size(3, control) != packed.size()) ITFAILS;

	int twodsize = packed.size() - sizeof(double);
	if (PT::get_packed_particle_size(2, control) != twodsize)      ITFAILS;
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

void transport_test()
{
    // make a six-cell mesh (defined by OS_Input)
    SP<rtt_imc_test::Parser> parser(new rtt_imc_test::Parser("OS_Input"));
    OS_Builder               mb(parser);
    SP<OS_Mesh>              mesh = mb.build_Mesh();

    // make an Opacity with sigma = 0.0
    SP<Opacity<OS_Mesh> > opacity;
    {
	OS_Mesh::CCSF<double> planck(mesh);
	opacity = new Opacity<OS_Mesh>(planck, planck, planck);
    }

    if (opacity->num_cells() != mesh->num_cells()) ITFAILS;

    // make a tally
    SP<Tally<OS_Mesh> > tally;

    // make a particle and transport it
    vector<double> r(2);
    vector<double> rf(2);
    vector<double> rfref(2);
    vector<double> o(3);
    Rnd_Control    control(seed);

    {
	// start in cell 1 and transport in the positive x direction until it
	// leaves the system 
	r[0] = -0.5;
	r[1] = -0.5;
	o[0] = 1.0;
	o[1] = 0.0;
	o[2] = 0.0;
	
	Sprng rnd = control.get_rn(1);

	// make particle
	PT p(r, o, 0.5, 1, rnd, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != PT::ESCAPE) ITFAILS;
	if (!soft_equiv(p.get_ew(), 0.5))     ITFAILS;
	
	// finishing position
	rfref[0] = 2.0;
	rfref[1] = -0.5;
	rf       = p.get_r();
	if (!soft_equiv(rf.begin(), rf.end(), rfref.begin(), rfref.end()))
	    ITFAILS;

	if (tally->get_accum_n_reflections() != 0) ITFAILS;
	if (tally->get_accum_n_escaped()     != 1) ITFAILS;
	if (tally->get_accum_n_bndcross()    != 3) ITFAILS;
    }

    {
	// start in cell 4 and transport in the positive x direction until it
	// hits the time-cutoff 1cm later
	r[0] = -0.5;
	r[1] = 2.0;
	o[0] = 1.0;
	o[1] = 0.0;
	o[2] = 0.0;
	
	Sprng rnd = control.get_rn(4);

	// make particle
	PT p(r, o, 0.5, 4, rnd, 1.0, 1.0);
	p.set_time_left(1.0/rtt_mc::global::c);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != PT::CENSUS) ITFAILS;
	
	// finishing position
	rfref[0] = 0.5;
	rfref[1] = 2.0;
	rf       = p.get_r();
	if (!soft_equiv(rf.begin(), rf.end(), rfref.begin(), rfref.end()))
	    ITFAILS;

	if (tally->get_accum_n_reflections() != 0) ITFAILS;
	if (tally->get_accum_n_escaped()     != 0) ITFAILS;
	if (tally->get_accum_n_bndcross()    != 1) ITFAILS;
	if (tally->get_new_ncen_tot()        != 1) ITFAILS;
    }

    {
	// start in cell 3 and transport off the low x boundary and back
	r[0] = 1.5;
	r[1] = -1.0/3.0;

	double tanphip  = -r[1]/2.5;
	double phi      = rtt_mc::global::pi - atan(tanphip);
	double cosphi   = cos(phi);
	double sinphi   = sqrt(1.0 - cosphi*cosphi);
	double costheta = 0.0;
	double sintheta = 1.0;
	o[0] = sintheta * cosphi;
	o[1] = sintheta * sinphi;
	o[2] = costheta;

	// check angle
	double norm = sqrt(o[0]*o[0] + o[1]*o[1] + o[2]*o[2]);
	if (!soft_equiv(norm, 1.0)) ITFAILS;
	
	Sprng rnd = control.get_rn(4);

	// make particle
	PT p(r, o, 0.5, 3, rnd, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != PT::ESCAPE) ITFAILS;
	
	// finishing position
	rfref[0] = 2.0;
	rfref[1] = 3.0 * tanphip;
	rf       = p.get_r();
	if (!soft_equiv(rf.begin(), rf.end(), rfref.begin(), rfref.end()))
	    ITFAILS;

	if (tally->get_accum_n_reflections() != 1) ITFAILS;
	if (tally->get_accum_n_escaped()     != 1) ITFAILS;
	if (tally->get_accum_n_bndcross()    != 6) ITFAILS;
	if (tally->get_new_ncen_tot()        != 0) ITFAILS;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Particle transport test ok.");
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
	transport_test();
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
