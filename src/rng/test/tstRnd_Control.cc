//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRnd_Control.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  9 13:59:31 2002
 * \brief  Rnd_Control test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "rng_test.hh"
#include "../Release.hh"
#include "../Random.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::soft_equiv;

using namespace std;

int seed = 2452423;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void control_test()
{
    // make a controller
    Rnd_Control control(seed);

    // checks
    if (control.get_number() != 1000000000) ITFAILS;
    if (control.get_seed()   != 2452423)    ITFAILS;
    if (control.get_num()    != 0)          ITFAILS;

    // make some random numbers
    Sprng r0  = control.get_rn();
    if (control.get_num()    != 1)          ITFAILS;
    Sprng r1  = control.get_rn();
    if (control.get_num()    != 2)          ITFAILS;
    Sprng r2  = control.get_rn();
    if (control.get_num()    != 3)          ITFAILS;

    Sprng rr2 = control.get_rn(2);
    if (control.get_num()    != 3)          ITFAILS;

    Sprng rr1 = control.get_rn(1);
    if (control.get_num()    != 2)          ITFAILS;

    control.set_num(0);
    Sprng rr0 = control.get_rn();
    if (control.get_num()    != 1)          ITFAILS;

    for (int i = 0; i < 100; i++)
    {
	double rn0  = r0.ran();
	double rrn0 = rr0.ran();
	double rn1  = r1.ran();
	double rrn1 = rr1.ran();
	double rn2  = r2.ran();
	double rrn2 = rr2.ran();

	if (!soft_equiv(rn0, rrn0))         ITFAILS;
	if (!soft_equiv(rn1, rrn1))         ITFAILS;
	if (!soft_equiv(rn2, rrn2))         ITFAILS;

	if (soft_equiv(rn0, rrn1))          ITFAILS;
	if (soft_equiv(rn1, rrn2))          ITFAILS;
	if (soft_equiv(rn2, rrn0))          ITFAILS;
    }
    
    vector<char> pack = r0.pack();
    if (control.get_size() != pack.size())  ITFAILS;

    if (rtt_rng_test::passed)
	PASSMSG("Rnd_Control simple test ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_rng::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	control_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstRnd_Control, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_rng_test::passed) 
    {
        cout << "**** tstRnd_Control Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstRnd_Control." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstRnd_Control.cc
//---------------------------------------------------------------------------//
