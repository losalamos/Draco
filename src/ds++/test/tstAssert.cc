//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstAssert.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 12:11:22 2003
 * \brief  Assertion tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Assert.hh"
#include "../Release.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// The way this test article works is that each of the DBC macros are tested
// in a seperate function.  A falst condition is asserted using each macro,
// and after this follows a throw.  Two catch clauses are available, one to
// catch an assertion object, and one to catch anything else.  By comparing
// the exception that is actually caught with the one that should be caught
// given the DBC setting in force, we can determine whether each test passes
// or fails.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Make sure we can differentiate betweeen a std::runtime_error and a
// rtt_dsxx::assertion. 
//---------------------------------------------------------------------------//

void t1()
{
    std::cout << "ta test: ";
    try 
    {
	throw std::runtime_error( "hello1" );
    } 
    catch( rtt_dsxx::assertion &a )
    {
	std::cout << "failed" << std::endl;
    }
    catch( ... )
    {
	std::cout << "passed" << std::endl;
    }
}

//---------------------------------------------------------------------------//
// Make sure we can catch a rtt_dsxx::assertion and extract the error
// message. 
//---------------------------------------------------------------------------//

void t2()
{
    std::cout << "t2-a test: ";
    std::string error_message;
    try 
    {
	throw rtt_dsxx::assertion( "hello1", "myfile", 42 );
    } 
    catch( rtt_dsxx::assertion &a )
    {
	std::cout << "passed" << std::endl;
	error_message = std::string( a.what() );
    }
    catch( ... )
    {
	std::cout << "failed" << std::endl;
    }

    // Make sure we can extract the error message.

    std::cout << "t2-b test: ";
    std::string const compare_value( 
	"Assertion: hello1, failed in myfile, line 42.\n" ); 
    if ( error_message.compare( compare_value ) == 0 )
	std::cout << "passed" << std::endl;
    else
	std::cout << "failed" << std::endl;
}

//---------------------------------------------------------------------------//
// Test throwing and catching of a literal
//---------------------------------------------------------------------------//

void t3()
{
    std::cout << "t3 test: ";
    try 
    {
	throw "hello";
    } 
    catch( rtt_dsxx::assertion &a )
    {
	std::cout << "failed" << std::endl;
    }
    catch( const char* msg )
    {
	std::cout << "passed" << std::endl;
    }
    catch( ... )
    {
	std::cout << "failed" << std::endl;
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Require() macro.
//---------------------------------------------------------------------------//

void trequire()
{
    std::cout << "t-Require test: ";
    try {
	Require( 0 );
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion& a )
    {
#if DBC & 1
	std::cout << "passed\n";
#else
	std::cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 1
	std::cout << "failed\n";
#else
	std::cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Check() macro.
//---------------------------------------------------------------------------//

void tcheck()
{
    std::cout << "t-Check test: ";
    try {
	Check( false );
	throw std::runtime_error( std::string( "tstAssert: t2()" ) );
    }
    catch( rtt_dsxx::assertion& a )
    {
#if DBC & 2
	std::cout << "passed\n";
#else
	std::cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 2
	std::cout << "failed\n";
#else
	std::cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Ensure() macro.
//---------------------------------------------------------------------------//

void tensure()
{
    int x = 0;
    Remember(x = 5);

    std::cout << "t-Ensure test: ";
    try {
	Ensure(0);
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion& a )
    {
#if DBC & 4
	std::cout << "passed\n";
#else
	std::cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 4
	std::cout << "failed\n";
#else
	std::cout << "passed\n";
#endif
    }

#if DBC & 4
    if (x != 5) ITFAILS;
#else
    if (x != 0) ITFAILS;
#endif
}

//---------------------------------------------------------------------------//
// Check the operation of the Assert() macro, which works like Check().
//---------------------------------------------------------------------------//

void tassert()
{
    std::cout << "t-Assert test: ";
    try {
	Assert(0);
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion& a )
    {
#if DBC & 2
	std::cout << "passed\n";
#else
	std::cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 2
	std::cout << "failed\n";
#else
	std::cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Basic test of the Insist() macro.
//---------------------------------------------------------------------------//

void tinsist()
{
    std::cout << "t-Insist test: ";
    try {
	Insist( 0, "You must be kidding!" );
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion& a ) {
	std::cout << "passed\n";
    }
    catch(...) {
	std::cout << "failed\n";
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}
    
    // >>> UNIT TESTS
    
    // Test basic throw and catch functionality.
    t1();
    t2();
    t3();

    // Test Design-by-Constract macros.
    trequire();
    tcheck();
    tensure();
    tassert();
    tinsist();

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstAssert Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstAssert." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstAssert.cc
//---------------------------------------------------------------------------//
