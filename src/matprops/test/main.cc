//----------------------------------*-C++-*----------------------------------//
// main.cc
// Randy M. Roberts
// Mon May  4 13:32:55 1998
//---------------------------------------------------------------------------//
// @> Driver for the matprops test routines.
//---------------------------------------------------------------------------//

#include "ds++/Assert.hh"
#include <iostream>
#include "testMmcMatProp.hh"
using std::cerr;
using std::endl;
using std::cout;
using std::string;

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main(int argc , char *argv[])
{
    try
    {

	for (int arg=1; arg < argc; arg++)
	    {
		if (std::string(argv[arg]) == "--version")
		    {
			version(argv[0]);
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
