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
    C4::Init( argc, argv );

    int node = C4::node();
    int nodes = C4::nodes();

    cout << "Hello from " << node << endl;

    C4::gsync();

    cout << "Hello again from " << node << endl;

    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tgsync.cc
//---------------------------------------------------------------------------//
