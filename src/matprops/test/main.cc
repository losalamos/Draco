//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Mon May  4 13:32:55 1998
//---------------------------------------------------------------------------//
// @> Driver for the matprops test routines.
//---------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include <iostream>
#include "matprops/test/testMmcMatProp.hh"
using std::cerr;
using std::endl;

int main()
{
    try
    {
    	void testBilinDoit();
	testBilinDoit();
	
	void testMatProp();
	testMatProp();
	
	testMmcMatProp xxx;
	xxx.execute_test();

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
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
