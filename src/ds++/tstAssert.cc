//----------------------------------*-C++-*----------------------------------//
// tstAssert.cc
// Geoffrey Furnish
// Mon Jul 28 13:26:55 1997
//---------------------------------------------------------------------------//
// @> Test code for checking behavior of the assertion facility.
//---------------------------------------------------------------------------//

#include "Assert.hh"
using namespace dsxx;

#include <iostream.h>

//---------------------------------------------------------------------------//
// Basic test of the assert macro.
//---------------------------------------------------------------------------//

void t1()
{
    try {
	Assert(0);
	cout << "t1 passed through \"Assert(0)\"." << endl;
    }
    catch( assertion& a ) {
	cout << "a.what() = " << a.what() << endl;
    }
    catch(...) {
	cout << "Big trouble in t1()!" << endl;
    }
}

//---------------------------------------------------------------------------//
// Basic test of the Insist macro.
//---------------------------------------------------------------------------//

void t2()
{
    try {
	Insist( 0, "You must be kidding!" );
	cout << "t2 passed through \"Insist(0,...)\"." << endl;
    }
    catch( assertion& a ) {
	cout << "a.what() = " << a.what() << endl;
    }
    catch(...) {
	cout << "Big trouble in t2()!" << endl;
    }
}

int main( int argc, char *argv[] )
{
    t1();
    t2();
}

//---------------------------------------------------------------------------//
//                              end of tstAssert.cc
//---------------------------------------------------------------------------//
