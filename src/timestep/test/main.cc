//----------------------------------*-C++-*----------------------------------//
// main.cc
// John McGhee
// Fri May  1 09:43:49 1998
//---------------------------------------------------------------------------//
// @> A driver for the time-step manager test facility.
//---------------------------------------------------------------------------//`
#include "timestep/ts_manager.hh"

#include "timestep/test/test_timestep.hh"

#include "ds++/Assert.hh"

#include "c4/global.hh"

#include <iostream>
using std::cerr;
using std::endl;

int main (int argc, char *argv[])
{
    try
    {
	C4::Init(argc, argv);

	// Create a nested scope since tester may be doing C4 stuff in its dtor.
	{
	    test_timestep tester;

	    tester.execute_test();
	}

	C4::Finalize();
	return 0;
    }
    catch (const dsxx::assertion &ass)
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
