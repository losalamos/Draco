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
#include "../Random_Walk.hh"
#include "../Release.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Constants.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_imc::Random_Walk;
using rtt_imc::Random_Walk_Sampling_Tables;
using rtt_dsxx::soft_equiv;
using rtt_dsxx::SP;

typedef rtt_mc::OS_Mesh MT;

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
// TESTS
//---------------------------------------------------------------------------//

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
    Random_Walk_Sampling_Tables table;
    {
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

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3)) ITFAILS;

	// a = 0.425, b = .3
	Ro = 0.1;
	t  = 0.01;
	D  = 0.425;
	R1 = 0.03;

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3)) ITFAILS;

	// a = 0.65, b= .9
	Ro = 0.1;
	t  = 0.01;
	D  = 0.65;
	R1 = 0.09;

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3)) ITFAILS;

	// a = 15.0, b = .1
	Ro = 0.1;
	t  = 0.01;
	D  = 15.0;
	R1 = 0.01;

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3)) ITFAILS;

	// a = 25.0, b = .2
	Ro = 0.1;
	t  = 0.01;
	D  = 25.0;
	R1 = 0.02;

	ran = get_prob_radius(t,D,Ro,R1) / get_prob_radius(t,D,Ro,Ro);
	R   = table.get_radius(t, D, Ro, ran); 

	if (!soft_equiv(R, R1, 1.e-3)) ITFAILS;
	
    }
    
    if (rtt_imc_test::passed)
	PASSMSG("Random Walk sampling tables for Rcen ok.");
}

//---------------------------------------------------------------------------//

void random_walk_OS_Mesh_test()
{
    Random_Walk<MT> rw;

    if (rtt_imc_test::passed)
	PASSMSG("Random Walk on OS_Mesh ok.");
}

//---------------------------------------------------------------------------//

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
	// >>> UNIT TESTS
	
	// sampling tables tests
	sampling_table_prob_exit_t_test();
	sampling_table_Rcen_test();

	// random walk tests
	random_walk_OS_Mesh_test();
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
