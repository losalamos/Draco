//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstMath.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:33:26 2001
 * \brief  Math functions test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_mc::global::soft_equiv;
using rtt_mc::global::integer_modulo;

//---------------------------------------------------------------------------//
// TESTS
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

    // test the integer modulo function
    if (integer_modulo(895,896) != 895)                      ITFAILS;
    if (integer_modulo(1999999999,2000000000) != 1999999999) ITFAILS;
    if (integer_modulo(2000000000,2000000000) != 0)          ITFAILS;
    if (integer_modulo(2000003091,2000000000) != 3091)       ITFAILS;
    
    // check macros
    if (INTEGER_MODULO_1E9(999999) != 999999) ITFAILS;
    if (INTEGER_MODULO_1E9(1000000000) != 0)  ITFAILS;
    if (INTEGER_MODULO_1E9(1000000002) != 2)  ITFAILS;
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

	// test math functions
	test_math();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstMath, " << ass.what()
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
	    cout << "**** tstMath Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstMath on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstMath.cc
//---------------------------------------------------------------------------//
