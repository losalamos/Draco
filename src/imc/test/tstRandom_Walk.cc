//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstRandom_Walk.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 28 13:32:27 2003
 * \brief  Random_Walk test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "imc_test.hh"
#include "IMC_Test.hh"
#include "../Random_Walk.hh"
#include "../Diffusion_Opacity.hh"
#include "../Fleck_Factors.hh"
#include "../Release.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Constants.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;

using rtt_imc_test::Parser;
using rtt_imc::Random_Walk;
using rtt_imc::Random_Walk_Sampling_Tables;
using rtt_imc::Fleck_Factors;
using rtt_imc::Diffusion_Opacity;
using rtt_mc::OS_Builder;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::SP;
using rtt_mc::global::pi;

typedef rtt_mc::OS_Mesh MT;
typedef pair<vector<double>, vector<double> > pair_dbl;

int seed = 234224;

//---------------------------------------------------------------------------//
// TESTING SERVICES
//---------------------------------------------------------------------------//

double integrate_prob_exit(double t, double D, double Ro)
{
    using rtt_mc::global::pi;

    double A = exp(-pi * pi / (Ro * Ro) * D * t);

    bool not_done = true;

    int n          = 1;
    double n2      = 0.0;
    double sum     = 0.0;
    double old_sum = 10.0;
    while (not_done)
    {
	n2   = static_cast<double>(n*n);
	sum += -2.0 * pow(A, n2) * pow(-1.0, static_cast<double>(n));

	if (soft_equiv(sum, old_sum, 1.0e-14)) not_done = false;

	old_sum = sum;
	
	n++;
    }

    // calc prob_exit
    sum = 1.0 - sum;

    Check (sum >= 0.0 && sum <= 1.0);
    return sum;
}

//---------------------------------------------------------------------------//

double get_prob_radius(double t, double D, double Ro, double R1)
{
    using rtt_mc::global::pi;

    double A = exp(-pi * pi / (Ro * Ro) * D * t);
    double B = R1/Ro;

    bool not_done = true;

    int n          = 1;
    double rn      = 1;
    double n2      = 0.0;
    double sum     = 0.0;
    double old_sum = 10.0;
    while (not_done)
    {
	rn   = static_cast<double>(n);
	n2   = rn*rn;
	sum += 2.0 * pow(A, n2) * (sin(rn*pi*B)/(rn*pi) - B * cos(rn*pi*B));

	if (soft_equiv(sum, old_sum, 1.0e-14)) not_done = false;

	old_sum = sum;
	
	n++;
    }

    Check (sum >= 0.0 && sum <= 1.0);
    return sum;
}

//---------------------------------------------------------------------------//

SP<Diffusion_Opacity<MT> > get_diff_opacity(SP<MT> mesh)
{
    // make data
    SP<Fleck_Factors<MT> > fleck(new Fleck_Factors<MT>(mesh));
    MT::CCSF<double> ross(mesh);

    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	fleck->fleck(i) = 1.0 / 2.5;
	ross(i)         = 100.0 + static_cast<double>(i);
    }

    // make the diffusion opacity
    SP<Diffusion_Opacity<MT> > diff(
	new Diffusion_Opacity<MT>(fleck, ross));
					
    return diff;
}

//===========================================================================//
// RANDOM WALK SAMPLING TABLE TESTS
//===========================================================================//

void sampling_table_prob_exit_t_test()
{
    Random_Walk_Sampling_Tables table;

    // test prob exit sampling
    {
	double Ro   = 0.0;
	double t    = 0.0;
	double D    = 0.0;
	double prob = 0.0;
	double ref  = 0.0;
	double t_e  = 0.0;

	// a = 0.01
	Ro = 0.1;
	t  = 0.01;
	D  = 0.01;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, 0.15671E-09)) ITFAILS;
	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.055
	Ro = 0.1;
	t  = 0.01;
	D  = 0.055;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, ref, 5.0e-2)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.1
	Ro = 0.1;
	t  = 0.01;
	D  = 0.1;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, 0.29290E+00)) ITFAILS;
	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.125
	Ro = 0.1;
	t  = 0.01;
	D  = 0.125;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.3
	Ro = 0.1;
	t  = 0.01;
	D  = 0.3;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, 0.89647E+00)) ITFAILS;
	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.425
	Ro = 0.1;
	t  = 0.01;
	D  = 0.425;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 0.65
	Ro = 0.1;
	t  = 0.01;
	D  = 0.65;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, t, 1.e-3)) ITFAILS;

	// a = 25.0
	Ro = 0.1;
	t  = 0.01;
	D  = 25.0;

	prob = table.get_prob_exit(t, D, Ro);
	ref  = integrate_prob_exit(t, D, Ro);

	if (!soft_equiv(prob, 1.0))         ITFAILS;
	if (!soft_equiv(prob, ref, 1.0e-3)) ITFAILS;

	t_e = table.get_elapsed_time(D, Ro, prob);

	if (!soft_equiv(t_e, 20.0 * Ro * Ro / D, 1.e-3)) ITFAILS;
    }

    if (rtt_imc_test::passed)
	PASSMSG("Random Walk sampling tables for prob_exit and t ok.");
}

//---------------------------------------------------------------------------//

void sampling_table_Rcen_test()
{
    cout << endl << "Entering R sampling test" << endl;

    Random_Walk_Sampling_Tables table;
    {
	cout << endl;

	double Ro  = 0;
	double t   = 0;
	double D   = 0;
	double R1  = 0;
	double ran = 0;
	double R   = 0;

	// a = .3, b = .7
	Ro = 0.1;
	t  = 0.01;
	D  = 0.3;
	R1 = 0.07;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = 0.425, b = .3
	Ro = 0.1;
	t  = 0.01;
	D  = 0.425;
	R1 = 0.03;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = 0.65, b= .9
	Ro = 0.1;
	t  = 0.01;
	D  = 0.65;
	R1 = 0.09;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran);

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	} 

	// a = 15.0, b = .1
	Ro = 0.1;
	t  = 0.01;
	D  = 15.0;
	R1 = 0.01;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = 25.0, b = .2
	Ro = 0.1;
	t  = 0.01;
	D  = 25.0;
	R1 = 0.02;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";
	
	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = 15.0, b = .15
	Ro = 0.1;
	t  = 0.01;
	D  = 15.0;
	R1 = 0.015;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = .65, b = .27
	Ro = 0.1;
	t  = 0.01;
	D  = 0.65;
	R1 = 0.027;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 5.e-2))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = .65, b = .62
	Ro = 0.1;
	t  = 0.01;
	D  = 0.65;
	R1 = 0.062;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	// a = .11, b = .21
	Ro = 0.1;
	t  = 0.01;
	D  = 0.11;
	R1 = 0.021;

	cout << "Checking for a=" << D*t/(Ro*Ro) << " and b="
	     << R1/Ro << " .... ";
	
	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 5.e-2))
	{
	    ITFAILS;
	}
	else 
	{
	    cout << "ok" << endl;
	}

	cout << endl;
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Random Walk sampling tables for Rcen ok.");
}

//---------------------------------------------------------------------------//

void interpolation_error_test()
{
    cout << endl;

    Random_Walk_Sampling_Tables table; 

    // test prob exit
    {
	double prob_exit = 0.0;
	double ref       = 0.0;
	double max_error = 0.0;
	double error     = 0.0;
	double D         = 0.1;
	double t         = 0.001;
	double Ro        = 0.0;
	double ma        = 0.0;
	
	double a = 0.04;
	while (a <= 20)
	{
	    // calc Ro
	    Ro = sqrt(D*t/a);
	    
	    // calc prob_exit
	    prob_exit = table.get_prob_exit(t, D, Ro);
	    ref       = integrate_prob_exit(t, D, Ro);
	    
	    // calculate difference
	    error     = fabs(prob_exit - ref) / ref;
	    max_error = max(error, max_error);
	    if (error == max_error) ma = a;
	    
	    // advance a
	    a += 0.04;
	}
	
	cout << "The maximum error interpolating Prob_exit is " << max_error
	     << " at a = " << ma << endl;
    }

    // test time of exit
    {
	double prob_exit = 0.0;
	double t         = 0.0;
	double max_error = 0.0;
	double error     = 0.0;
	double D         = 0.1;
	double ref       = 0.001;
	double Ro        = 0.001;
	double mpe       = 0.0;

	prob_exit = 0.05;
	while (prob_exit <= 1.0)
	{
	    // get a value for t
	    t = table.get_elapsed_time(D, Ro, prob_exit);

	    // now integrate to get prob_exit
	    ref = integrate_prob_exit(t, D, Ro);

	    // estimate an error
	    error     = fabs(prob_exit - ref) / ref;
	    max_error = max(error, max_error);
	    if (error == max_error) mpe = prob_exit;

	    // advance prob_exit
	    prob_exit += 0.05;
	}

	cout << "The maximum error interpolating t exit is " << max_error
	     << " at Prob_exit = " << mpe << endl;
    }

    // test of R1
    {
	double a = 0.04;
	double b = 0.025;

	double t   = 0.001;
	double Ro  = 0.0;
	double D   = 0.035;
	double ran = 0.0;
	double R1  = 0.0;
	double ref = 0.0;

	double error     = 0.0;
	double max_error = 0.0;
	double ma        = 0.0;
	double mb        = 0.0;

	while (a <= 20.0)
	{
	    b = 0.025;
	    
	    // calculate Ro
	    Ro = sqrt(t*D/a);

	    while (b <= 1.0)
	    {
		// calc R1
		ref = b * Ro;

		// calculate ran
		ran = get_prob_radius(t,D,Ro,ref) / 
		    get_prob_radius(t,D,Ro,Ro);
		
		// calculate R1
		R1 = table.get_radius(t,D,Ro,ran);

		// estimate an error
		error     = fabs(R1 - ref) / ref;
		max_error = max(error, max_error);
		if (error == max_error) 
		{
		    ma = a;
		    mb = b;
		}

		// advance b
		b += 0.025;
	    }

	    // advance a 
	    a += 0.04; 
	}

	cout << "The maximum error interpolating R1 is " << max_error
	     << " at a = " << ma << " and b = " << mb << endl << endl;
    }
}

//===========================================================================//
// RANDOM WALK TESTS
//===========================================================================//

void random_walk_test()
{
    // build a 3D mesh
    SP<MT> mesh;
    {
	SP<Parser> p(new Parser("OS_Input_3D"));
	OS_Builder b(p);
	mesh = b.build_Mesh();

	if (mesh->num_cells() != 12) ITFAILS;
    }
    
    // get the diffusion opacity
    SP<Diffusion_Opacity<MT> > diff = get_diff_opacity(mesh);

    // make random walk
    Random_Walk<MT> rw(mesh, diff);

    // make two random number generators
    Rnd_Control control(seed);
    Sprng rng  = control.get_rn(11);
    Sprng rngr = control.get_rn(11);

    // random walk table
    Random_Walk_Sampling_Tables table;

    // particle data
    vector<double> r(3);
    vector<double> o(3);

    // check do_a_random_walk function
    // Ross = 100 + cell_index
    {
	cout << "Random Walk example in cell 1" << endl;
	cout << "=============================" << endl;

	double radius = 0.0;
	double d_col  = 0.0;
	double d_cen  = 0.0;
	double telap  = 0.0;
	double tleft  = 0.0;
	double D      = 0.0;
	bool   toc    = false;
       
	double P_exit = 0.0;
	double ran    = 0.0;
	double d      = 0.0;
	double dref   = 0.0;
	double cost   = 0.0;
	double phi    = 0.0;
	vector<double> oldr(3);
	vector<double> refo(3);
	
	// >>> cell 1
	radius = 1.0 / 10.0;
	d_col  = 1.0 / 102.0;
	d_cen  = 1.0 / 5.0;

	if (!rw.do_a_random_walk(1, radius, d_col, d_cen)) ITFAILS;
	
	// now do the random walk
	r[0]  = -0.5;
	r[1]  = 0.1;
	r[2]  = 0.6;
	tleft = 0.001;

	oldr = r;

	// check it
	D      = diff->get_random_walk_D(1);
	P_exit = table.get_prob_exit(tleft, D, radius);
	ran    = rngr.ran();

	cout << "Probability of escape is : " << P_exit << endl;
	cout << "Random number is         : " << ran << endl;

	// the particle will go to census
	d = table.get_radius(tleft, D, radius, rngr.ran());
	
	// do the random walk
	telap = rw.random_walk(r, o, tleft, 1, rng, toc);

	if (!toc)                      ITFAILS;
	if (!soft_equiv(tleft, 0.0))   ITFAILS;
	if (!soft_equiv(telap, 0.001)) ITFAILS;

	// check particle position
	dref = 0.0;
	for (int i = 0; i < 3; i++)
	    dref += (r[i] - oldr[i]) * (r[i] - oldr[i]);
	dref = sqrt(dref);
	if (!soft_equiv(d, dref)) ITFAILS;

	pair_dbl rn = mesh->sample_pos_on_sphere(1, oldr, d, rngr);
	oldr = rn.first;
	if (!soft_equiv(r.begin(), r.end(), oldr.begin(), oldr.end())) ITFAILS;

	// check particle direction
	refo = rn.second;
	cost = sqrt(rngr.ran());
	phi  = 2.0 * pi * rngr.ran();
	mesh->get_Coord().calc_omega(cost, phi, refo);

	if (!soft_equiv(o.begin(), o.end(), refo.begin(), refo.end())) ITFAILS;

	if (!soft_equiv(rng.ran(), rngr.ran()))                        ITFAILS;

	cout << "=============================" << endl << endl;

	// cell 1 where random walk will be off
	radius = 1.0 / 10.0;
	d_col  = 1.0 / 9.9;
	d_cen  = 1.0 / 5.0;

	if (rw.do_a_random_walk(1, radius, d_col, d_cen)) ITFAILS;

	radius = 1.0 / 10.0;
	d_col  = 1.0 / 102.0;
	d_cen  = 1.0 / 11.0;

	if (rw.do_a_random_walk(1, radius, d_col, d_cen)) ITFAILS;

	radius = 1.0 / 101.5;
	d_col  = 1.0 / 102.0;
	d_cen  = 1.0 / 10.0;

	if (rw.do_a_random_walk(1, radius, d_col, d_cen)) ITFAILS;

	// >>> cell 11
	cout << "Random Walk example in cell 11" << endl;
	cout << "==============================" << endl;

	radius = 1.0 / 10.0;
	d_col  = 1.0 / 20.0;
	d_cen  = 1.0 / 5.0;

	if (!rw.do_a_random_walk(11, radius, d_col, d_cen)) ITFAILS;
	
	// now do the random walk
	r[0]  = 0.5;
	r[1]  = 2.1;
	r[2]  = 1.77;
	tleft = 0.001;

	oldr = r;

	// check it
	D      = diff->get_random_walk_D(11);
	P_exit = table.get_prob_exit(tleft, D, radius);
	ran    = rngr.ran();

	cout << "Probability of escape is : " << P_exit << endl;
	cout << "Random number is         : " << ran << endl;

	// the particle makes it to the surface of the sphere
	telap = rw.random_walk(r, o, tleft, 11, rng, toc);

	if (toc)                               ITFAILS;
	if (!soft_equiv(telap + tleft, 0.001)) ITFAILS;

	// check time elapsed
	double rtelap = table.get_elapsed_time(D, radius, ran);
	if (!soft_equiv(telap, rtelap)) ITFAILS;
	if (telap < 0.0)                ITFAILS;
	if (tleft < 0.0)                ITFAILS;

	// check particle position
	dref = 0.0;
	for (int i = 0; i < 3; i++)
	    dref += (r[i] - oldr[i]) * (r[i] - oldr[i]);
	dref = sqrt(dref);
	if (!soft_equiv(radius, dref)) ITFAILS;

	rn   = mesh->sample_pos_on_sphere(11, oldr, radius, rngr);
	oldr = rn.first;
	if (!soft_equiv(r.begin(), r.end(), oldr.begin(), oldr.end())) ITFAILS;

	// check particle direction
	refo = rn.second;
	cost = sqrt(rngr.ran());
	phi  = 2.0 * pi * rngr.ran();
	mesh->get_Coord().calc_omega(cost, phi, refo);

	if (!soft_equiv(o.begin(), o.end(), refo.begin(), refo.end())) ITFAILS;

	if (!soft_equiv(rng.ran(), rngr.ran()))                        ITFAILS;

	cout << "==============================" << endl << endl;

	// try some where it is off
	radius = 1.0 / 10.0;
	d_col  = 1.0 / 9.9;
	d_cen  = 1.0 / 5.0;

	if (rw.do_a_random_walk(11, radius, d_col, d_cen)) ITFAILS;

	radius = 1.0 / 10.0;
	d_col  = 1.0 / 102.0;
	d_cen  = 1.0 / 11.0;

	if (rw.do_a_random_walk(11, radius, d_col, d_cen)) ITFAILS;

	radius = 1.0 / 111.5;
	d_col  = 1.0 / 112.0;
	d_cen  = 1.0 / 10.0;

	if (rw.do_a_random_walk(11, radius, d_col, d_cen)) ITFAILS;
    }
   
    if (rtt_imc_test::passed)
	PASSMSG("Random Walk on OS_Mesh ok.");
}

//===========================================================================//
// MAIN
//===========================================================================//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_imc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {

	// only run on 1 processor
	if (rtt_c4::node() == 0)
	{
	    // >>> UNIT TESTS
	
	    // sampling tables tests
	    sampling_table_prob_exit_t_test();
	    sampling_table_Rcen_test();
	    interpolation_error_test();

	    // random walk tests
	    random_walk_test();
	}
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstRandom_Walk, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstRandom_Walk Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstRandom_Walk on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstRandom_Walk.cc
//---------------------------------------------------------------------------//
