//----------------------------------*-C++-*----------------------------------//
// tstNestMap.cc
// Shawn Pautz
// Tue Mar 30 08:22:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../NestMap.t.hh"
#include <iostream>

#ifdef SDP

using rtt_dsxx::NestMap;
using std::cout;
using std::endl;

bool passed = true;

template class NestMap<std::string, int>;

#endif

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    std::cout << progname << ": version " << version << std::endl;
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
    
#ifdef SDP

    cout << "Initiating test of the NestMap class.\n";

    NestMap<std::string, int> intTree;

    intTree.open("one", 1);   // stack: 1
    intTree.close();          // stack:
    intTree.open("one", 1);   // stack: 1
    intTree.close();          // stack:
    intTree.open("two", 2);   // stack: 2
    intTree.close();          // stack:

    intTree.open("one", 1);            // stack: 1
    intTree.open("one one", 11);       // stack: 1 11
    intTree.open("one one one", 111);  // stack: 1 11 111
    intTree.close();                   // stack: 1 11
    intTree.close();                   // stack: 1
    intTree.open("one two", 12);       // stack: 1 12
    intTree.close();                   // stack: 1
    intTree.open("one three", 13);     // stack: 1 13
    intTree.close();                   // stack: 1
    intTree.close();                   // stack:
    intTree.open("three", 3);          // stack: 3
    intTree.open("three three", 33);   // stack: 3 33
    intTree.close();                   // stack: 3
    intTree.close();                   // stack:

    cout << "\nDumping objects" << endl;
    for (NestMap<std::string, int>::iterator iter = intTree.begin();
	 iter != intTree.end(); ++iter)
    {
	cout << (*iter).first << endl;
    }

// Print the status of the test.

    cout << endl;
    cout <<     "***********************************" << endl;
    if (passed) 
    {
        cout << "**** NestMap Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** NestMap Self Test: FAILED ****" << endl;
    }
    cout <<     "***********************************" << endl;
    cout << endl;

    cout << "Done testing NestMap class.\n";

#endif

    std::cout << "**** NestMap Self Test: PASSED ****" << std::endl;

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstNestMap.cc
//---------------------------------------------------------------------------//
