//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Mon May  4 13:32:55 1998
//---------------------------------------------------------------------------//
// @> Driver for the matprops test routines.
//---------------------------------------------------------------------------//

#include "c4/global.hh"

#include "ds++/Assert.hh"
#include <iostream>
#include "testMmcMatProp.hh"
#include "../Release.hh"

using std::cerr;
using std::endl;
using std::cout;
using std::string;

void version(const std::string &progname)
{
    cout << progname << ": version " << rtt_matprops::release() << endl;
}

int main(int argc , char *argv[])
{
    try
    {
	C4::Init( argc, argv );
	
	for (int arg=1; arg < argc; arg++)
	{
	    if (std::string(argv[arg]) == "--version")
	    {
		if (C4::node() == 0)
		    version(argv[0]);
		C4::Finalize();
		return 0;
	    }
	}

    	void testBilinDoit();
	testBilinDoit();
	
	void testMatProp();
	testMatProp();
	
	testMmcMatProp xxx;
	xxx.execute_test();

    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cerr << "Test: FAILED: assertion: " << ass.what() << endl;
    }
    catch (const std::exception &ass)
    {
	cerr << "Test: FAILED: exception: " << ass.what() << endl;
    }
    catch (...)
    {
	cerr << "Test: FAILED: unknown exception" << endl;
    }

    cout << "Done testing matprops.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
