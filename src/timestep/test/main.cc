//----------------------------------*-C++-*----------------------------------//
// main.cc
// John McGhee
// Fri May  1 09:43:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//`
#include "timestep/ts_manager.hh"

#include "timestep/test/test_timestep.hh"

#include "ds++/Assert.hh"

#include <iostream>
using std::cerr;
using std::endl;

int main ()
{
    try
    {
	test_timestep tester;

	tester.execute_test();

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
