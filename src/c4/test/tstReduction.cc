//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstReduction.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 15:41:00 2002
 * \brief  C4 Reduction test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include "../Release.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_c4::global_sum;
using rtt_c4::global_prod;
using rtt_c4::global_min;
using rtt_c4::global_max;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void elemental_reduction()
{
    // test ints
    int xint = rtt_c4::node() + 1;
    global_sum(xint);

    int int_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
	int_answer += i + 1;

    if (xint != int_answer) ITFAILS;

    // test longs
    long xlong = rtt_c4::node() + 1000;
    global_sum(xlong);

    long long_answer = 0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
	long_answer += i + 1000;

    if (xlong != long_answer) ITFAILS;

    // test doubles
    double xdbl = static_cast<double>(rtt_c4::node()) + 0.1;
    global_sum(xdbl);

    double dbl_answer = 0.0;
    for (int i = 0; i < rtt_c4::nodes(); i++)
	dbl_answer += static_cast<double>(i) + 0.1;

    if (!soft_equiv(xdbl, dbl_answer)) ITFAILS;

    // test product
    xlong = rtt_c4::node() + 1;
    global_prod(xlong);

    long_answer = 1;
    for (int i = 0; i < rtt_c4::nodes(); i++)
	long_answer *= (i + 1);

    if (xlong != long_answer) ITFAILS;

    // test min
    xdbl = 0.5 + rtt_c4::node();
    global_min(xdbl);
    
    if (!soft_equiv(xdbl, 0.5)) ITFAILS;

    // test max
    xdbl = 0.7 + rtt_c4::node();
    global_max(xdbl);

    if (!soft_equiv(xdbl, rtt_c4::nodes() - 0.3)) ITFAILS;

    if (rtt_c4_test::passed)
	PASSMSG("Elemental reductions ok.");
}

//---------------------------------------------------------------------------//

void array_reduction()
{
    // make a vector of doubles
    vector<double> x(100);
    vector<double> prod(100, 1.0);
    vector<double> sum(100, 0.0);
    vector<double> min(100, 0.0);
    vector<double> max(100, 0.0);
    
    // fill it
    for (int i = 0; i < 100; i++)
    {
	x[i]  = rtt_c4::node() + 0.11;
	for (int j = 0; j < rtt_c4::nodes(); j++)
	{
	    sum[i]  += (j + 0.11);
	    prod[i] *= (j + 0.11);
	}
	min[i] = 0.11;
	max[i] = rtt_c4::nodes() + 0.11 - 1.0;
    }

    vector<double> c;

    c = x;
    global_sum(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), sum.begin(), sum.end())) ITFAILS;

    c = x;
    global_prod(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), prod.begin(), prod.end())) ITFAILS;

    c = x;
    global_min(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), min.begin(), min.end())) ITFAILS;

    c = x;
    global_max(&c[0], 100);
    if (!soft_equiv(c.begin(), c.end(), max.begin(), max.end())) ITFAILS;
    
    if (rtt_c4_test::passed)
	PASSMSG("Array reductions ok.");
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
		cout << argv[0] << ": version " << rtt_c4::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
       
	elemental_reduction();
	array_reduction();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstReduction, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_c4_test::passed) 
	{
	    cout << "**** tstReduction Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstReduction on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstReduction.cc
//---------------------------------------------------------------------------//
