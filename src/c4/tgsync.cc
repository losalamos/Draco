//----------------------------------*-C++-*----------------------------------//
// tgsync.cc
// Geoffrey Furnish
// Wed Jan 18 13:39:39 1995
//---------------------------------------------------------------------------//
// @> Test program for C4_gsync().
//---------------------------------------------------------------------------//

#include "c4/global.hh"

#include <iostream.h>

//---------------------------------------------------------------------------//
// Verify that global synchronization is working.  Two groups of messages are
// printed, they should not overlap, but ordering within each group is
// nondeterministic. 
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    C4_Init( argc, argv );

    int node = C4_node();
    int nodes = C4_nodes();

    cout << "Hello from " << node << endl;

    C4_gsync();

    cout << "Hello again from " << node << endl;

    C4_Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tgsync.cc
//---------------------------------------------------------------------------//
