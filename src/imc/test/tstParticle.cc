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
#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "../Opacity.hh"
#include "../Tally.hh"
#include "../Global.hh"
#include "../Frequency.hh"
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

using rtt_imc::Gray_Particle;
using rtt_imc::Multigroup_Particle;
using rtt_imc::Opacity;
using rtt_imc::Tally;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_rng::Sprng;
using rtt_rng::Rnd_Control;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::SP;

// random number seed
int seed = 395731;

typedef Gray_Particle<OS_Mesh>        GPT;
typedef Multigroup_Particle<OS_Mesh>  MGPT;
typedef rtt_imc::Gray_Frequency       G;
typedef rtt_imc::Multigroup_Frequency MG;

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
    vector<double> ref4(100);
    {
	Sprng r1 = control.get_rn(1);
	Sprng r2 = control.get_rn(2);
	Sprng r3 = control.get_rn(3);
	Sprng r4 = control.get_rn(4);

	for (int i = 0; i < 100; i++)
	{
	    ref1[i] = r1.ran();
	    ref2[i] = r2.ran();
	    ref3[i] = r3.ran();
	    ref4[i] = r4.ran();
	}
    }

    // make gray and multigroup particles
    vector<double> r(3);
    vector<double> o(3);
    double         ew   = 0.45;
    double         tl   = 0.15;
    double         fr   = 0.5;
    int            cell = 10;
    int            grp  = 2;
    Sprng          rnd1 = control.get_rn(1);
    Sprng          rnd2 = control.get_rn(2);
	
    // throw some random numbers
    for (int i = 0; i < 15; i++)
    {
	rnd1.ran();
	rnd2.ran();
    }

    r[0] = 0.5;
    r[1] = 0.75;
    r[2] = 4.0;
    o[0] = 0.5;
    o[1] = 0.5;
    o[2] = std::sqrt(2.0)/2.0;
    GPT  gray(r, o, ew, cell, rnd1, fr, tl);
    MGPT mg(r, o, ew, cell, rnd2, grp, fr, tl);

    // output particles
    if (C4::node() == 0)
    {
	cout << gray << endl;

	cout << mg << endl;
    }

    // check the particle after a copy
    {
	GPT gp   = gray;
	MGPT mgp = mg;

	if (gp != gray) ITFAILS;
	if (mgp != mg)  ITFAILS;

	// check it
	vector<double> rf(3);
	vector<double> of(3);

	rf[0] = 0.5;
	rf[1] = 0.75;
	rf[2] = 4.0;
	of[0] = 0.5;
	of[1] = 0.5;
	of[2] = std::sqrt(2.0)/2.0;
	
	vector<double> gr  = gp.get_r();
	vector<double> go  = gp.get_omega();
	vector<double> mgr = mgp.get_r();
	vector<double> mgo = mgp.get_omega();
	
	if (!soft_equiv(gr.begin(), gr.end(), rf.begin(), rf.end()))   ITFAILS;
	if (!soft_equiv(go.begin(), go.end(), of.begin(), of.end()))   ITFAILS;
	if (!soft_equiv(mgr.begin(), mgr.end(), rf.begin(), rf.end())) ITFAILS;
	if (!soft_equiv(mgo.begin(), mgo.end(), of.begin(), of.end())) ITFAILS;
	if (gp.get_cell() != 10)                                       ITFAILS;
	if (!soft_equiv(gp.get_ew(), 0.45))                            ITFAILS;
	if (!gp.status())                                              ITFAILS;
	if (gp.get_descriptor() != GPT::BORN)                          ITFAILS;
	if (mgp.get_cell() != 10)                                      ITFAILS;
	if (!soft_equiv(mgp.get_ew(), 0.45))                           ITFAILS;
	if (!mgp.status())                                             ITFAILS;
	if (mgp.get_descriptor() != MGPT::BORN)                        ITFAILS;
	if (mgp.get_group_index() != 2)                                ITFAILS;

	// check random numbers; the stream is the same as the original
	// particle 
	for (int i = 0; i < 20; i++)
	{
	    double r   = gp.get_random().ran();
	    double ref = ref1[i+15];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
	for (int i = 0; i < 20; i++)
	{
	    double r   = mgp.get_random().ran();
	    double ref = ref2[i+15];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }

    // check that the original particles are on the same stream
    for (int i = 0; i < 65; i++)
    {
	double r   = gray.get_random().ran();
	double ref = ref1[i+35];
	if (!soft_equiv(r, ref)) ITFAILS;
    }
    for (int i = 0; i < 65; i++)
    {
	double r   = mg.get_random().ran();
	double ref = ref2[i+35];
	if (!soft_equiv(r, ref)) ITFAILS;
    }

    // check set functions
    {
	Sprng r3 = control.get_rn(3);
	Sprng r4 = control.get_rn(4);

	gray.set_random(r3);
	mg.set_random(r4);

	for (int i = 0; i < 50; i++)
	{
	    double r   = gray.get_random().ran();
	    double ref = ref3[i];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}

	for (int i = 0; i < 50; i++)
	{
	    double r   = mg.get_random().ran();
	    double ref = ref4[i];
	    if (!soft_equiv(r, ref)) ITFAILS;
	}
    }

    // gray checks
    {
	gray.set_cell(12);
	gray.set_ew(0.25);

	if (gray.get_cell() != 12)               ITFAILS;
	if (!soft_equiv(gray.get_ew(), 0.25))    ITFAILS;

	gray.kill_particle();

	if (gray.status())                       ITFAILS;

	gray.reset_status();

	if (!gray.status())                      ITFAILS;

	if (gray.get_descriptor() != GPT::BORN)  ITFAILS;
	gray.set_descriptor(GPT::CENSUS);
	if (gray.get_descriptor() != GPT::CENSUS)ITFAILS;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Gray_Particle simple tests ok.");

    // multigroup checks
    {
	mg.set_cell(12);
	mg.set_ew(0.25);

	if (mg.get_cell() != 12)                 ITFAILS;
	if (!soft_equiv(mg.get_ew(), 0.25))      ITFAILS;

	mg.kill_particle();

	if (mg.status())                         ITFAILS;

	mg.reset_status();

	if (!mg.status())                        ITFAILS;

	if (mg.get_descriptor() != MGPT::BORN)   ITFAILS;
	mg.set_descriptor(MGPT::ESCAPE);
	if (mg.get_descriptor() != MGPT::ESCAPE) ITFAILS;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Multigroup_Particle simple tests ok.");
}

//---------------------------------------------------------------------------//

void gray_pack_test()
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
	GPT p(r, o, ew, cell, rnd, fr, tl);

	// pack the particle
	packed = p.pack();

	// throw more random numbers
	rnd.ran();

	// check the size packed state function
	if (GPT::get_packed_particle_size(2, control) == packed.size()) ITFAILS;
	if (GPT::get_packed_particle_size(3, control) != packed.size()) ITFAILS;

	int twodsize = packed.size() - sizeof(double);
	if (GPT::get_packed_particle_size(2, control) != twodsize)      ITFAILS;
    }

    // unpack it and test
    {
	GPT p(packed);

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
	if (p.get_descriptor() != GPT::UNPACKED)                   ITFAILS;

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
	PASSMSG("Gray_Particle pack test ok.");
}

//---------------------------------------------------------------------------//

void mg_pack_test()
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
	int            grp  = 3;
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
	MGPT p(r, o, ew, cell, rnd, grp, fr, tl);

	// pack the particle
	packed = p.pack();

	// throw more random numbers
	rnd.ran();

	// check the size packed state function
	if (MGPT::get_packed_particle_size(2, control) == packed.size()) ITFAILS;
	if (MGPT::get_packed_particle_size(3, control) != packed.size()) ITFAILS;

	int twodsize = packed.size() - sizeof(double);
	if (MGPT::get_packed_particle_size(2, control) != twodsize)      ITFAILS;
	
	int diff = MGPT::get_packed_particle_size(2, control) -
	    GPT::get_packed_particle_size(2, control);
	if (diff != sizeof(int))                                         ITFAILS;
    }

    // unpack it and test
    {
	MGPT p(packed);

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
	if (p.get_group_index() != 3)                              ITFAILS;
	if (p.get_descriptor() != GPT::UNPACKED)                   ITFAILS;

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
	PASSMSG("Multigroup_Particle pack test ok.");
}

//---------------------------------------------------------------------------//

void gray_transport_test()
{
    // make a six-cell mesh (defined by OS_Input)
    SP<rtt_imc_test::Parser> parser(new rtt_imc_test::Parser("OS_Input"));
    OS_Builder               mb(parser);
    SP<OS_Mesh>              mesh = mb.build_Mesh();

    // make an Opacity with sigma = 0.0
    SP<Opacity<OS_Mesh, G> > opacity;
    {
	SP<G> frequency(new G);
	OS_Mesh::CCSF<double> planck(mesh);
	OS_Mesh::CCSF<double> fleck(mesh);
	for (int i = 1; i <= mesh->num_cells(); i++)
	    fleck(i) = 1.0;

	opacity = new Opacity<OS_Mesh, G>(frequency,planck, planck, fleck);
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
	GPT p(r, o, 0.5, 1, rnd, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                        ITFAILS;
	if (p.get_descriptor() != GPT::ESCAPE) ITFAILS;
	if (!soft_equiv(p.get_ew(), 0.5))      ITFAILS;
	
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
	GPT p(r, o, 0.5, 4, rnd, 1.0, 1.0);
	p.set_time_left(1.0/rtt_mc::global::c);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != GPT::CENSUS) ITFAILS;
	
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
	GPT p(r, o, 0.5, 3, rnd, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != GPT::ESCAPE) ITFAILS;
	
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
	PASSMSG("Gray_Particle transport test ok.");
}

//---------------------------------------------------------------------------//

void mg_transport_test()
{
    // make a six-cell mesh (defined by OS_Input)
    SP<rtt_imc_test::Parser> parser(new rtt_imc_test::Parser("OS_Input"));
    OS_Builder               mb(parser);
    SP<OS_Mesh>              mesh = mb.build_Mesh();

    // make an Opacity with sigma = 0.0
    SP<Opacity<OS_Mesh, MG> > opacity;
    {
	vector<double> grps(3);
	grps[0] = 0.01;
	grps[1] = 0.1;
	grps[2] = 1.0;
	SP<MG> frequency(new MG(grps));

	// make cross sections
	OS_Mesh::CCSF<vector<double> > planck(mesh);
	OS_Mesh::CCSF<double>          int_planck(mesh);
	vector<double>                 xs(2, 0.0);

	for (int i = 1; i <= mesh->num_cells(); i++)
	{
	    int_planck(i) = 1.0;
	    planck(i)     = xs;
	}

	opacity = new Opacity<OS_Mesh, MG>(frequency, planck, planck,
					   int_planck, int_planck, planck);
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
	MGPT p(r, o, 0.5, 1, rnd, 2, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != GPT::ESCAPE) ITFAILS;
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
	MGPT p(r, o, 0.5, 4, rnd, 1, 1.0, 1.0);
	p.set_time_left(1.0/rtt_mc::global::c);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != GPT::CENSUS) ITFAILS;
	
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
	MGPT p(r, o, 0.5, 3, rnd, 1, 1.0, 1.0);

	// make tally
	tally = new Tally<OS_Mesh>(mesh);

	p.transport(*mesh, *opacity, *tally);

	// check the particle
	if (p.status())                       ITFAILS;
	if (p.get_descriptor() != GPT::ESCAPE) ITFAILS;
	
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
	PASSMSG("Multigroup_Particle transport test ok.");
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

	// gray tests
	gray_pack_test();
	gray_transport_test();

	// multigroup tests
	mg_pack_test();
	mg_transport_test();
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
