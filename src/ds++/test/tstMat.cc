//----------------------------------*-C++-*----------------------------------//
// tstMat.cc
// Geoffrey Furnish
// Wed Apr  2 12:48:48 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::flush;

// This is heinous beyond description, but I /really/ want to test all the
// parts of the Matrix classes here.  Wantom abuse of CPP like that which
// follows, is just to prove that I mean business!

#define private public
#define protected public

#include "../Mat.hh"

using namespace rtt_dsxx;

template<class T, class A>
void ikv1( Mat1<T,A>& x )
{
    cout << "In ikv1." << endl;
}

void t1()
{
    cout << "t1: beginning.\n";

    {
    // Test Mat1<int>;

	Mat1<int> x;
	x.redim(5);

	ikv1( x );
    }

    cout << "t1 test: passed\n";
}

//---------------------------------------------------------------------------//
// Test the Mat2, every way we can think of.
//---------------------------------------------------------------------------//

void t2()
{
    cout << "t2: beginning.\n";
    
    {
	cout << "Testing default ctor. ";

	Mat2<int> x;

	cout  << "Done." << endl;
    }
    cout << "t2a test: passed.\n";

    {
	cout << "Testing conventional ctor. " << flush;;

    // Check various fundamental computations.
	Mat2<int> x(3,3);

	Assert( x.nx() == 3 );
	Assert( x.ny() == 3 );
	Assert( x.index(0,0) == 0 );
	Assert( x.size() == 9 );

	int k=0;
	for( int j=0; j < 3; j++ )
	    for( int i=0; i < 3; i++ )
		x(i,j) = k++;

	k = 0;
	for( Mat2<int>::iterator xi = x.begin(); xi != x.end(); )
	    if (*xi++ != k++)
		throw "bogus";

	cout << "Done." << endl;
    }
    cout << "t2b test: passed.\n";

    {
	cout << "Testing ctor with Bounds. " << flush;;

	Mat2<int> x( Bounds(1,3), Bounds(1,3) );

	Assert( x.nx() == 3 );
	Assert( x.ny() == 3 );
        Assert( x.index(1,1) == 4 );
	Assert( x.size() == 9 );

	int k=0;
	for( int j=1; j < 4; j++ )
	    for( int i=1; i < 4; i++ )
		x(i,j) = k++;

	k = 0;
	for( Mat2<int>::iterator xi = x.begin(); xi != x.end(); )
	    if (*xi++ != k++)
		throw "bogus";

	cout << "Done." << endl;
    }
    cout << "t2c test: passed.\n";

    cout << "Done testing Mat2<T>." << endl;
}


void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{

    try {
	for (int arg=1; arg < argc; arg++)
	    {
		if (std::string(argv[arg]) == "--version")
		    {
			version(argv[0]);
			return 0;
		    }
	    }
	
	cout << "Initiating test of the Mat family.\n";

	t1();
	t2();
    }
    catch( assertion& a )
    {
	cout << "Test: Failed assertion: " << a.what() << endl;
    }

    cout << "Done testing Mat family.\n";
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstMat.cc
//---------------------------------------------------------------------------//
