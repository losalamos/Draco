//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstShared_Lib.cc
 * \author Rob Lowrie
 * \date   Thu Apr 15 23:03:32 2004
 * \brief  Tests Shared_Lib
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

#include "../Shared_Lib.hh"
#include "../Assert.hh"
#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "ds_test.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Only compile if we have dlopen support on this platform
#ifdef NO_DLOPEN

void test_simple()
{
    PASSMSG("dlopen not configured/available for this platform.");
}

#else

void test_simple()
{
    // This lib must contain pow(double, double)
    string math_lib_name("/usr/lib/libm.so");

    Shared_Lib math(math_lib_name);
    
    if ( math.get_file_name() != math_lib_name ) ITFAILS;
    if ( not math.is_open() ) ITFAILS;
    
    // FP is a pointer to a function that returns a double and takes two
    // doubles as arguments.
    typedef double (* FP)(double, double);

    // Grab pow(double, double)
    FP my_pow = math.get_function<FP>("pow");

    const double x = 2.34342;

    if ( not soft_equiv(my_pow(x, 2), x * x) ) ITFAILS;

    { // Test copy ctor
	Shared_Lib sc(math);

	if ( sc.get_file_name() != math_lib_name ) ITFAILS;
	FP p = sc.get_function<FP>("pow");
	if ( not soft_equiv(p(x, 2), x * x) ) ITFAILS;
    }

    { // Test assignment
	Shared_Lib sc;

	sc = math;

	if ( sc.get_file_name() != math_lib_name ) ITFAILS;
	FP p = sc.get_function<FP>("pow");
	if ( not soft_equiv(p(x, 2), x * x) ) ITFAILS;
    }

    // Make sure everything is still working with original copy
    
    if ( math.get_file_name() != math_lib_name ) ITFAILS;
    if ( not math.is_open() ) ITFAILS;
    FP another_pow = math.get_function<FP>("pow");
    if ( not soft_equiv(another_pow(x, 2), x * x) ) ITFAILS;

    { // check get_function of a function not in the library
        bool caught = false;
        try
        {
	    FP no = math.get_function<FP>("not_in_math");
        }
        catch ( const rtt_dsxx::assertion &ass )
        {
            caught = true;
            std::ostringstream m;
            m << "Excellent! Caught assertion for get_function().";
            PASSMSG(m.str());
        }

        if ( not caught ) ITFAILS;
    }

#ifdef REQUIRE_ON
    { // check open() of an empty file name
        bool caught = false;
        try
        {
	    Shared_Lib t;
	    t.open("");
        }
        catch ( const rtt_dsxx::assertion &ass )
        {
            caught = true;
            std::ostringstream m;
            m << "Excellent! Caught assertion for open().";
            PASSMSG(m.str());
        }

        if ( not caught ) ITFAILS;
    }
#endif

    // Done testing

    if ( rtt_ds_test::passed )
    {
        PASSMSG("test_simple() ok.");
    }
}

#endif // NO_DLOPEN

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_dsxx::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	test_simple();
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tstShared_Lib, " 
		  << err.what()
		  << std::endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstShared_Lib, " 
		  << "An unknown exception was thrown."
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_ds_test::passed) 
    {
        std::cout << "**** tstShared_Lib Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstShared_Lib." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstShared_Lib.cc
//---------------------------------------------------------------------------//
