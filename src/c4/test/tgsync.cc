//----------------------------------*-C++-*----------------------------------//
// tgsync.cc
// Geoffrey Furnish
// Wed Jan 18 13:39:39 1995
//---------------------------------------------------------------------------//
// @> Test program for C4_gsync().
//---------------------------------------------------------------------------//

#include "../global.hh"

#include <iostream>
#include <string>

using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
// Verify that global synchronization is working.  Two groups of messages are
// printed, they should not overlap, but ordering within each group is
// nondeterministic. 
//---------------------------------------------------------------------------//

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

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

    C4::gsync();

    cout << "Hello again from " << node << endl;

    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tgsync.cc
//---------------------------------------------------------------------------//
