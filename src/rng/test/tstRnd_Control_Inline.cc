//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRnd_Control_Inline.cc
 * \author Paul Henning
 * \brief  Rnd_Control test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "rng_test.hh"
#include "../Release.hh"
#include "../Random_Inline.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>


using rtt_rng::Rnd_Control;
using rtt_rng::LF_Gen;


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
    LF_Gen r0; control.initialize(r0);
    if (control.get_num()    != 1)          ITFAILS;
    LF_Gen r1; control.initialize(r1);
    if (control.get_num()    != 2)          ITFAILS;
    LF_Gen r2; control.initialize(r2);
    if (control.get_num()    != 3)          ITFAILS;

    LF_Gen rr2; control.initialize(2, rr2);
    if (control.get_num()    != 3)          ITFAILS;

    LF_Gen rr1; control.initialize(1, rr1);
    if (control.get_num()    != 2)          ITFAILS;

    control.set_num(0);

    LF_Gen rr0; control.initialize(rr0);
    if (control.get_num()    != 1)          ITFAILS;

    for (int i = 0; i < 100; i++)
    {
	double rn0  = r0.ran();
	double rrn0 = rr0.ran();
	double rn1  = r1.ran();
	double rrn1 = rr1.ran();
	double rn2  = r2.ran();
	double rrn2 = rr2.ran();

	if (rn0 != rrn0)         ITFAILS;
	if (rn1 != rrn1)         ITFAILS;
	if (rn2 != rrn2)         ITFAILS;

	if (rn0 == rrn1)          ITFAILS;
	if (rn1 == rrn2)          ITFAILS;
	if (rn2 == rrn0)          ITFAILS;
    }
    

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
