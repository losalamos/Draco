//----------------------------------*-C++-*----------------------------------//
// tBSwap.cc
// Geoffrey Furnish
// Wed Jan 18 13:39:39 1995
//---------------------------------------------------------------------------//
// @> Test program for BSwap<T>.
//---------------------------------------------------------------------------//

#include "../global.hh"
#include "../BSwap.hh"
using namespace C4;

#include <iostream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
// This program demonstrates the use of the BSwap<T> class by passing an
// integer across all the nodes, doing a little work on it each time.  This
// is a very trivial example, but should be sufficient to test the class.
//---------------------------------------------------------------------------//

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}


int main( int argc, char *argv[] )
{
    Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
	{
	    if (std::string(argv[arg]) == "--version")
		{
		    version(argv[0]);
		    C4::Finalize();
		    return 0;
		}
	}

    int node = C4::node();
    int nodes = C4::nodes();

    cout << "Hello from " << node << endl;

    int i;

    if (node) {
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
	i=2;
	if (nodes > 1) {
	    BSwap<int> right( 1 );

	    right.send(i);
	} else
	    cout << "2 to the 1 power is " << i << endl;
    }

    if (node == nodes-1) {
    // We're the last node, so its our job to assess success/failure.

        if ( i == pow(2.,nodes) )
            cout << nodes << " Node BSwap Test: passed\n";
        else
            cout << nodes << " Node BSwap Test: failed\n";
    }

    C4::gsync();
    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tBSwap.cc
//---------------------------------------------------------------------------//
