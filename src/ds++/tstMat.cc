//----------------------------------*-C++-*----------------------------------//
// tstMat.cc
// Geoffrey Furnish
// Wed Apr  2 12:48:48 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "Mat.hh"

template<class T, class A>
void ikv1( Mat1<T,A>& x )
{
    cout << "In ikv1." << endl;
}

void t1()
{
    {
    // Test Mat1<int>;

	Mat1<int> x;
	x.redim(5);

	ikv1( x );
    }
}

void main( int argc, char *argv[] )
{
    cout << "Initiating test of the Mat family.\n";

    t1();

    cout << "Done testing Mat family.\n";
}

//---------------------------------------------------------------------------//
//                              end of tstMat.cc
//---------------------------------------------------------------------------//
