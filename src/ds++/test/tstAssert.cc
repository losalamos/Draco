//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstAssert.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 12:11:22 2003
 * \brief  Assertion tests.
 * \note   Copyright (c) 1997-2010 Los Alamos National Security, LLC
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

static void t1()
{
    std::cout << "t1 test: ";
    try 
    {
	throw std::runtime_error( "hello1" );
    } 
    catch( rtt_dsxx::assertion const & a )
    {
	std::cout << a.what() << std::endl;
	std::cout << "failed" << std::endl;
    }
    catch( ... )
    {
	std::cout << "passed" << std::endl;
    }
    return;
}

//---------------------------------------------------------------------------//
// Make sure we can catch a rtt_dsxx::assertion and extract the error
// message. 
//---------------------------------------------------------------------------//

static void t2()
{
    std::cout << "t2-a test: ";
    std::string error_message;
    try 
    {
	throw rtt_dsxx::assertion( "hello1", "myfile", 42 );
    } 
    catch( rtt_dsxx::assertion const & a )
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
    return;
}

//---------------------------------------------------------------------------//
// Test throwing and catching of a literal
//
//lint -e1775  do not warn about "const char*" not being a declared exception
//             type.  
//lint -e1752  do not warn about catching "const char*" instead of
//             catching a reference to an exception (see "More Effective C++"
//             for details about catching exceptions)
//---------------------------------------------------------------------------//

static void t3()
{
    std::cout << "t3 test: ";
    try 
    {
	throw "hello";
    } 
    catch( rtt_dsxx::assertion const & a )
    {
	std::cout << a.what() << std::endl;
	std::cout << "failed" << std::endl;
    }
    catch( const char* msg )
    {
	std::cout << "passed   "
		  << "msg = " << msg << std::endl;
    }
    catch( ... )
    {
	std::cout << "failed" << std::endl;
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the toss_cookies function.
// This function builds an error message and throws an exception.
//---------------------------------------------------------------------------//
static void ttoss_cookies()
{
    {
        std::cout << "ttoss_cookies test: ";
        try
        {
            std::string const msg("testing toss_cookies()");
            std::string const file("DummyFile.ext");
            int const line( 55 );
            rtt_dsxx::toss_cookies( msg, file, line );
            throw "Bogus!";
        }
        catch( rtt_dsxx::assertion const & /* error */ )
        {
            std::cout << "passed" << std::endl;
        }
        catch( ... )
        {
            std::cout << "failed" << std::endl;
        }
    }
    {
        std::cout << "ttoss_cookies_ptr test: ";
        try
        {
            char const * const msg("testing toss_cookies_ptr()");
            char const * const file("DummyFile.ext");
            int const line( 56 );
            rtt_dsxx::toss_cookies_ptr( msg, file, line );
            throw "Bogus!";
        }
        catch( rtt_dsxx::assertion const &  /* error */ )
        {
            std::cout << "passed" << std::endl;
        }
        catch( ... )
        {
            std::cout << "failed" << std::endl;
        } 
    }
    return;
}



//---------------------------------------------------------------------------//
// Check the operation of the Require() macro.
//---------------------------------------------------------------------------//

static void trequire()
{
    std::cout << "t-Require test: ";
    try {
	Require( 0 );
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion const & a )
    {
#if DBC & 1
	std::cout << "passed" << std::endl;

	std::cout << "t-Require message value test: ";
	{
	    std::string msg( a.what() );
	    std::string expected_value( "Assertion: 0, failed in" );
	    string::size_type idx = msg.find( expected_value );
	    if( idx != string::npos )
	    {
		cout << "passed" << std::endl;
	    }
	    else
	    {
		cout << "failed" << std::endl;
	    }
	}
	
#else
	std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
#endif
    }
    catch(...)
    {
#if DBC & 1
	std::cout << "failed" << std::endl;
#else
	std::cout << "passed" << std::endl;
#endif
    }
    return;
}

//---------------------------------------------------------------------------//
// Check the operation of the Check() macro.
//---------------------------------------------------------------------------//

static void tcheck()
{
    std::cout << "t-Check test: ";
    try {
	Check( false );
	throw std::runtime_error( std::string( "tstAssert: t2()" ) );
    }
    catch( rtt_dsxx::assertion const & a )
    {
#if DBC & 2
	std::cout << "passed" << std::endl;

	std::cout << "t-Check message value test: ";
	{
	    std::string msg( a.what() );
	    std::string expected_value( "Assertion: false, failed in" );
	    string::size_type idx = msg.find( expected_value );
	    if( idx != string::npos )
	    {
		cout << "passed" << std::endl;
	    }
	    else
	    {
		cout << "failed" << std::endl;
	    }
	}
#else
	std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
	std::string msg( a.what() );
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
    return;
}

//---------------------------------------------------------------------------//
// Check the operation of the Ensure() macro.
//---------------------------------------------------------------------------//

static void tensure()
{
    std::cout << "t-Ensure test: ";
    try {
	Ensure(0);
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion const & a )
    {
#if DBC & 4
	std::cout << "passed" << std::endl;

	std::cout << "t-Ensure message value test: ";
	{
	    std::string msg( a.what() );
	    std::string expected_value( "Assertion: 0, failed in" );
	    string::size_type idx = msg.find( expected_value );
	    if( idx != string::npos )
	    {
		cout << "passed" << std::endl;
	    }
	    else
	    {
		cout << "failed" << std::endl;
	    }
	}

#else
	std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
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
    return;
}

static void tremember()
{
    //lint -e774  do not warn about if tests always evaluating to False.  The
    //            #if confuses flexelint here.

    std::cout << "t-Remember test: ";

    int x = 0;
    Remember(x = 5);
#if DBC & 4
    if (x != 5) 
	std::cout << "failed" << std::endl;
    else
	std::cout << "passed" << std::endl;
#else
    if (x != 0) 
	std::cout << "failed" << std::endl;
    else
	std::cout << "passed" << std::endl;
#endif
    return;}

//---------------------------------------------------------------------------//
// Check the operation of the Assert() macro, which works like Check().
//---------------------------------------------------------------------------//

static void tassert()
{
    std::cout << "t-Assert test: ";
    try {
	Assert(0);
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion const & a )
    {
#if DBC & 2
	std::cout << "passed" << std::endl;

	std::cout << "t-Assert message value test: ";
	{
	    std::string msg( a.what() );
	    std::string expected_value( "Assertion: 0, failed in" );
	    string::size_type idx = msg.find( expected_value );
	    if( idx != string::npos )
	    {
		cout << "passed" << std::endl;
	    }
	    else
	    {
		cout << "failed" << std::endl;
	    }
	}
#else
	std::cout << "failed" << "\t" << "a.what() = " << a.what() << std::endl;
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
    return;}

//---------------------------------------------------------------------------//
// Basic test of the Insist() macro.
//---------------------------------------------------------------------------//

static void tinsist()
{
    //lint -e506  Do not warn about constant value boolean in the Insist
    //            test. 
    {
        std::cout << "t-Insist test: ";
        std::string insist_message( "You must be kidding!" );
        try {
            Insist( 0, insist_message );
            throw "Bogus!";
        }
        catch( rtt_dsxx::assertion const & a ) 
        {
            std::cout << "passed" << std::endl;
            
            std::cout << "t-Insist message value test: ";
            {
                bool passed( true );
                std::string msg( a.what() );
                std::string expected_value( "You must be kidding!" );
                string::size_type idx( msg.find( expected_value ) );
                if( idx == string::npos ) passed=false;
                idx = msg.find( insist_message );
                if( idx == string::npos ) passed=false;
                if( passed )
                    cout << "passed" << std::endl;
                else
                    cout << "failed" << std::endl;
            }
        }
        catch(...) 
        {
            std::cout << "failed" << std::endl;
        }
    }
    
    {
        std::cout << "t-Insist ptr test: ";
        char const * const insist_message( "You must be kidding!" );
        try {
            Insist_ptr( 0, insist_message );
            throw "Bogus!";
        }
        catch( rtt_dsxx::assertion const & a ) 
        {
            std::cout << "passed" << std::endl;
            
            std::cout << "t-Insist ptr message value test: ";
            {
                bool passed( true );
                std::string msg( a.what() );
                std::string expected_value( "You must be kidding!" );
                string::size_type idx( msg.find( expected_value ) );
                if( idx == string::npos ) passed=false;
                idx = msg.find( insist_message );
                if( idx == string::npos ) passed=false;
                if( passed )
                    cout << "passed" << std::endl;
                else
                    cout << "failed" << std::endl;
            }
        }
        catch(...) 
        {
            std::cout << "failed" << std::endl;
        }
    }
    
    return;
}

//---------------------------------------------------------------------------//
// Basic test of the Insist_ptr() macro.
//---------------------------------------------------------------------------//

static void tinsist_ptr()
{
    //lint -e506  Do not warn about constant value boolean in the Insist
    //            test. 
    
    std::cout << "t-Insist test: ";
    try {
	Insist( 0, "You must be kidding!" );
	throw "Bogus!";
    }
    catch( rtt_dsxx::assertion const & a ) 
    {
	std::cout << "passed" << std::endl;
        
	std::cout << "t-Insist_ptr message value test: ";
	{
	    bool passed( true );
	    std::string msg( a.what() );
	    std::string expected_value( "You must be kidding!" );
	    string::size_type idx( msg.find( expected_value ) );
	    if( idx == string::npos ) passed=false;
	    if( passed )
		cout << "passed" << std::endl;
	    else
		cout << "failed" << std::endl;
	}
    }
    catch(...) 
    {
	std::cout << "failed" << std::endl;
    }
    return;
}

void tverbose_error()
{
    std::string const message( rtt_dsxx::verbose_error(
                             std::string("This is an error.") ) );
    std::cout << "verbose_error() test: ";
    if( message.find( std::string("Host")) != std::string::npos &&
        message.find( std::string("PID") ) != std::string::npos )
    {
        cout << "passed" << std::endl;
        rtt_ds_test::passed = rtt_ds_test::passed && true;
    }
    else
    {
        cout << "failed" << std::endl;
        rtt_ds_test::passed = rtt_ds_test::passed && false;
    }
    
    return;
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    //lint -e30 -e85 -e24 -e715 -e818 Suppress warnings about use of argv 
    //          (string comparison, unknown length, etc.)

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if( string( argv[arg] ).find( "--version" ) == 0 )
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

    // Test mechanics of Assert funtions.
    ttoss_cookies();
    
    // Test Design-by-Constract macros.
    trequire();
    tcheck();
    tensure();
    tremember();
    tassert();
    tinsist();
    tinsist_ptr();

    // fancy ouput
    tverbose_error();
    
    // status of test
    cout <<     "\n*********************************************\n";
    if (rtt_ds_test::passed) 
        cout << "**** tstAssert Test: PASSED\n";
    cout <<     "*********************************************\n\n"
         << "Done testing tstAssert." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstAssert.cc
//---------------------------------------------------------------------------//
