//----------------------------------*-C++-*----------------------------------//
// tBSwap.cc
// Geoffrey Furnish
// Wed Jan 18 13:39:39 1995
//---------------------------------------------------------------------------//
// @> Test program for BSwap<T>.
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/NodeInfo.hh"
#include "c4/BSwap.hh"

#include <iostream.h>

#ifdef __GNUC__
#include "c4/BSwap.cc"
template class BSwap<int>;
#endif

#ifdef _POWER
#include "c4/BSwap.cc"
#pragma define ( BSwap<int> )
#endif

#ifdef __DECCXX
#include "c4/BSwap.cc"
#pragma define_template BSwap<int>
#endif

using namespace C4;

//---------------------------------------------------------------------------//
// This program demonstrates the use of the BSwap<T> class by passing an
// integer across all the nodes, doing a little work on it each time.  This
// is a very trivial example, but should be sufficient to test the class.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    Init( argc, argv );

    int node = C4::node();
    int nodes = C4::nodes();

    cout << "Hello from " << node << endl;

    if (node) {
	int i;
	BSwap<int> left( node-1 );
	left.recv(i);
	i *= 2;
	if (node < nodes-1) {
	    BSwap<int> right( node+1 );
	    right.send(i);
	} else {
	    cout << "2 to the " << nodes << " power is " << i << endl;
	}

    } else {
	int i=2;
	if (nodes > 1) {
	    BSwap<int> right( 1 );

	    right.send(i);
	} else
	    cout << "2 to the 1 power is " << i << endl;
    }

    C4::gsync();
    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tBSwap.cc
//---------------------------------------------------------------------------//
