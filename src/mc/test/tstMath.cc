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
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_mc::global::integer_modulo;
using rtt_mc::global::linear_interpolate;
using rtt_mc::global::log_log_interpolate;
using rtt_dsxx::soft_equiv;

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

void test_interpolate()
{
    // function y = 2.5 * x - 1.0

    // define boundary points
    double x1 = 1.0;
    double y1 = 2.5 * x1 - 1.0;
    double x2 = 3.0;
    double y2 = 2.5 * x2 - 1.0;

    double x   = 1.452;
    double y   = linear_interpolate(x1, x2, y1, y2, x);
    double ref = 2.5 * x - 1.0;

    if (!soft_equiv(y, ref)) ITFAILS;

    // try another one
    x1 = 1.45;
    y1 = 2.5 * x1 - 1.0;
    x2 = 1.1;
    y2 = 2.5 * x2 - 1.0;

    x   = 1.33;
    y   = linear_interpolate(x1, x2, y1, y2, x);
    ref = 2.5 * x - 1.0;

    if (!soft_equiv(y, ref)) ITFAILS;
 
    if (rtt_mc_test::passed)
	PASSMSG("Linear interpolation checks ok.");

    // function y = 1.1 * x^(-3.2)

    x1 = 1.0;
    y1 = 1.1 * pow(x1, -3.2);
    x2 = 1.2;
    y2 = 1.1 * pow(x2, -3.2);

    x   = 1.11;
    y   = log_log_interpolate(x1, x2, y1, y2, x);
    ref = 1.1 * pow(x, -3.2);

    if (!soft_equiv(y, ref)) ITFAILS;

    x1 = 3.31;
    y1 = 1.1 * pow(x1, -3.2);
    x2 = 2.12;
    y2 = 1.1 * pow(x2, -3.2);

    x   = 2.22;
    y   = log_log_interpolate(x1, x2, y1, y2, x);
    ref = 1.1 * pow(x, -3.2);

    if (!soft_equiv(y, ref)) ITFAILS;
 
    if (rtt_mc_test::passed)
	PASSMSG("Log-log interpolation checks ok.");
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

	// interpolate test
	test_interpolate();
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
