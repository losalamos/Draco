//----------------------------------*-C++-*----------------------------------//
// tstAssert.cc
// Geoffrey Furnish
// Mon Jul 28 13:26:55 1997
//---------------------------------------------------------------------------//
// @> Test code for checking behavior of the assertion facility.
//---------------------------------------------------------------------------//

#include "../Assert.hh"
using namespace rtt_dsxx;

#include <iostream>
#include <string>

using std::cout;
using std::endl;

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
// Check the operation of the Require() macro.
//---------------------------------------------------------------------------//

void t1()
{
    cout << "t1 test: ";
    try {
	Require(0);
	throw "Bogus!";
    }
    catch( assertion& a )
    {
#if DBC & 1
	cout << "passed\n";
#else
	cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 1
	cout << "failed\n";
#else
	cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Check() macro.
//---------------------------------------------------------------------------//

void t2()
{
    cout << "t2 test: ";
    try {
	Check(0);
	throw "Bogus!";
    }
    catch( assertion& a )
    {
#if DBC & 2
	cout << "passed\n";
#else
	cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 2
	cout << "failed\n";
#else
	cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Ensure() macro.
//---------------------------------------------------------------------------//

void t3()
{
    cout << "t3 test: ";
    try {
	Ensure(0);
	throw "Bogus!";
    }
    catch( assertion& a )
    {
#if DBC & 4
	cout << "passed\n";
#else
	cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 4
	cout << "failed\n";
#else
	cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Check the operation of the Assert() macro, which works like Check().
//---------------------------------------------------------------------------//

void t4()
{
    cout << "t4 test: ";
    try {
	Assert(0);
	throw "Bogus!";
    }
    catch( assertion& a )
    {
#if DBC & 2
	cout << "passed\n";
#else
	cout << "failed\n";
#endif
    }
    catch(...)
    {
#if DBC & 2
	cout << "failed\n";
#else
	cout << "passed\n";
#endif
    }
}

//---------------------------------------------------------------------------//
// Basic test of the Insist() macro.
//---------------------------------------------------------------------------//

void t5()
{
    cout << "t5 test: ";
    try {
	Insist( 0, "You must be kidding!" );
	throw "Bogus!";
    }
    catch( assertion& a ) {
	cout << "passed\n";
    }
    catch(...) {
	cout << "failed\n";
    }
}

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{

    for (int arg=1; arg < argc; arg++)
	{
	    if (std::string(argv[arg]) == "--version")
		{
		    version(argv[0]);
		    return 0;
		}
	}

    t1();
    t2();
    t3();
    t4();
    t5();
}

//---------------------------------------------------------------------------//
//                              end of tstAssert.cc
//---------------------------------------------------------------------------//
