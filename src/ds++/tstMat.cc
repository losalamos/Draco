//----------------------------------*-C++-*----------------------------------//
// tstMat.cc
// Geoffrey Furnish
// Wed Apr  2 12:48:48 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include <iostream>

using namespace std;

// This is heinous beyond description, but I /really/ want to test all the
// parts of the Matrix classes here.  Wantom abuse of CPP like that which
// follows, is just to prove that I mean business!

#define private public
#define protected public

#include "Mat.hh"

using namespace dsxx;

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

    cout << "t1: passed\n";
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
    cout << "t2a: passed.\n";

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
    cout << "t2b: passed.\n";

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
    cout << "t2c: passed.\n";

    cout << "Done testing Mat2<T>." << endl;
}

int main( int argc, char *argv[] )
{
    cout << "Initiating test of the Mat family.\n";

    try {
	t1();
	t2();
    }
    catch( assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

    cout << "Done testing Mat family.\n";
}

//---------------------------------------------------------------------------//
//                              end of tstMat.cc
//---------------------------------------------------------------------------//
