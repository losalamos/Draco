//----------------------------------*-C++-*----------------------------------//
/*!
 * \file fourier_test.cc
 * \author John Gulick
 * \date  Wed Jan 5 13:01:10 2000
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iomanip>
#include "../fourier.hh"

bool passed = true;

void version(const std::string &progname)
{
  std::string version = "1.0.0";
  cout << progname << ": version " << version << endl;
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
    cout << "Initiating test of Fourier Analysis package.\n";
    Fourier test_fourier;
    test_fourier.input();
    test_fourier.solve();
    test_fourier.print();
    cout << endl;
    cout <<   "******************************************" << endl;
    if (passed) 
    {
	cout << "**** Fourier Analysis Package Self Test: PASSED ****" << endl;
    }
    else
    {
	cout << "**** Fourier Analysis Package Self Test: FAILED ****" << endl;
    }
    cout <<   "******************************************" << endl;
    cout << endl;
    cout << "Done testing Fourier Analysis package.\n";
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of fourier_test.cc
//---------------------------------------------------------------------------//
