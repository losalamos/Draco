//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstMath.cc
 * \author Thomas M. Evans
 * \date   Thu Jan 13 14:51:33 2000
 * \brief  Math Functions test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "c4/global.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

using rtt_mc::global::soft_equiv;
using rtt_mc::global::mod_with_2e9;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void test_math()
{
    double x = 0.9876543212345678;
    double y = 0.9876543212345678;

    if (!soft_equiv(x, y, 1.e-16)) ITFAILS;
    if (!soft_equiv(x, y))         ITFAILS;

    double z = 0.9876543212345679;

    if (soft_equiv(x, z, 1.e-16)) ITFAILS;

    double a = 0.987654321234;
    
    if (!soft_equiv(x, a)) ITFAILS;

    a = 0.987654321233;

    if (soft_equiv(x, a)) ITFAILS;      

    // checks for the new "reference=zero" coding 4aug00
    double zero = 0.0;
    if ( soft_equiv( 1.0e-10, zero)) ITFAILS;
    if ( soft_equiv(-1.0e-10, zero)) ITFAILS;
    if (!soft_equiv(-1.0e-35, zero)) ITFAILS;
    if (!soft_equiv( 1.0e-35, zero)) ITFAILS;

    // test the mod_with_2e9 function
    if (mod_with_2e9(895)        != 895)        ITFAILS;
    if (mod_with_2e9(1999999999) != 1999999999) ITFAILS;
    if (mod_with_2e9(2000000000) != 0)          ITFAILS;
    if (mod_with_2e9(2000003091) != 3091)       ITFAILS;

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
	    cout << argv[0] << ": version " << rtt_mc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// test math functions
	test_math();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "You are a bona-fide Uncle-#@&%er: " << ass.what() << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "********************************" << endl;
    if (passed) 
    {
        cout << "**** Math Self Test: PASSED ****" << endl;
    }
    cout <<     "********************************" << endl;
    cout << endl;

    cout << "Done testing Math." << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstMath.cc
//---------------------------------------------------------------------------//
