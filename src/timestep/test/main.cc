//----------------------------------*-C++-*----------------------------------//
// main.cc
// John McGhee
// Fri May  1 09:43:49 1998
//---------------------------------------------------------------------------//
// @> A driver for the time-step manager test facility.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ts_manager.hh"
#include "../Release.hh"

#include "test_timestep.hh"

#include "ds++/Assert.hh"

#include "c4/global.hh"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::string;

void version(const std::string &progname)
{
    cout << progname << ": version " << rtt_timestep::release() << endl;
}

int main (int argc, char *argv[])
{

    try
    {
	C4::Init(argc, argv);

	for (int arg=1; arg < argc; arg++)
	    {
		if (std::string(argv[arg]) == "--version")
		    {
			version(argv[0]);
			C4::Finalize();
			return 0;
		    }
	    }

	// Create a nested scope since tester may be doing C4 stuff in its dtor.
	{
	    test_timestep tester;

	    tester.execute_test();
	}

	C4::Finalize();
	return 0;
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cerr << "assert failed: " << ass.what() << endl;
    }
    catch (const std::exception &ass)
    {
	cerr << "exception: " << ass.what() << endl;
    }
    catch (...)
    {
	cerr << "unknown exception" << endl;
    }
    return 1;
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
