//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Fri Mar 20 12:04:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testP13T.hh"
#include "ds++/Assert.hh"

#include <stdexcept>
#include <iostream>
using std::cerr;
using std::endl;

int main()
{
    try
    {
	cerr << "In main()" << endl;

	XTM::testP13T test;

	cerr << "before test.solve()" << endl;
    
	test.solve();

	cerr << "after test.solve()" << endl;
    }
    catch (const char *str)
    {
	cerr << "caught: " << str << endl;
	return 1;
    }
    catch (const dsxx::assertion &ass)
    {
	cerr << "caught assertion exception: " << ass.what() << endl;
	return 1;
    }
    catch (const std::runtime_error &rtx)
    {
	cerr << "caught assertion exception: " << rtx.what() << endl;
	return 1;
    }
    catch (...)
    {
	cerr << "caught unknown exception" << endl;
	return 1;
    }
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
