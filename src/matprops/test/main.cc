//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Mon May  4 13:32:55 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include <iostream>
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
