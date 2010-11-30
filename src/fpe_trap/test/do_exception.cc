//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/test/do_exception.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:33:59 2005
 * \brief  Does a floating-point exception.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../fpe_trap.hh"
#include "../Release.hh"
#include <ds++/Assert.hh>
#include <fstream>
#include <cmath>
#include <iostream>

int main(int argc, char *argv[]);

using namespace std;

/*
  Usage: do_exception test
  
  If test is 0, then simple floating point operations
  are done which should not cause an error.
     
  Otherwise, other test values should cause an exception.
  Specifically, valid test values are
     1: test double division by zero
     2: test sqrt(-1)
     3: test overflow

  The file output.dat documents what happened during all tests.
*/

//---------------------------------------------------------------------------//

int run_test(int /*argc*/, char *argv[])
{
    std::ofstream f;
    f.open("output.dat");

    if ( rtt_fpe_trap::enable_fpe() )
    {
        // Platform supported.
        f << "- fpe_trap: This platform is supported" << endl;
    }
    else
    {
        // Platform not supported.
        f << "- fpe_trap: This platform is not supported\n";
        f.close();
        return 0;
    }

    // Accept a command line argument with value 1, 2 or 3.
    int test(-101);
    sscanf(argv[1], "%d", &test);
    Insist(test >= 0 && test <= 3, "Bad test value.");

    double zero(    0.0 ); // for double division by zero
    double neg(    -1.0 ); // for sqrt(-1.0)
    double result( -1.0 );

    // Certain tests may be optimized away by the compiler, by recogonizing
    // the constants set above and precomputing the results below.  So do
    // something here to hopefully avoid this.  This tricks the optimizer, at
    // least for gnu and KCC.

    if ( test < -100 )
    { // this should never happen
	Insist(0, "Something is very wrong.");
	zero = neg = 1.0; // trick the optimizer?
    }

    switch ( test )
    {
        case 0:
            // The test_filter.py triggers on the keyword 'signal', so I will
            // use '5ignal' instead.
            f << "- Case zero: this operation should not throw a 5ignal."
              << " The result should be 2." << endl;
            result = 1.0 + zero + sqrt(-neg);
            f << "  result = " << result << endl;
            break;
        case 1:
            f << "- Trying a div_by_zero operation" << endl;
            result = 1.0 / zero; // should fail here
            f << "  result = " << result << endl;
            break;
        case 2:
            f << "- Trying to evaluate sqrt(-1.0)" << endl;
            result = sqrt(neg); // should fail here
            f << "  result = " << result << endl;
            break;
        case 3:
        {
            f << "- Trying to cause an overflow condition" << endl;
            result = 2.0;
            for ( size_t i = 0; i < 100; i++ )
            {
                // exp() should work, but on 64-bit linux, it's not raising the
                // overflow flag:
                // result = exp(result); // should fail at some i
                
                // ... so instead:
                result = result * result * result * result
                         * result; // should fail at some i
                
            }
            f << "  result = " << result << endl;
            break;
        }
    }
    // close the log file.
    f.close();
    return 0;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    Insist(argc == 2, "Wrong number of args.");

    // banner
    cout << "This is " << argv[0] << ": version "
         << rtt_fpe_trap::release() << endl; 
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	    return 0;
    
    try
    {
        run_test(argc,argv);
    }
    catch (exception &err)
    {
        if ( rtt_fpe_trap::enable_fpe() )
        {
            // keyword 'signal' shows up as a failure when processed by
            // test_filter.py 
            cout << "While running " << argv[0] << ", " 
                 << "a 5ignal was successfully caught.\n\t"
                // << err.what()
                 << endl;
            return 0;
        }
        else
        {
            cout << "ERROR: While running " << argv[0] << ", "
                 << "An exception was caught when it was not expected.\n\t"
                // << err.what()
                 << endl;
        }        
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown."
             << endl;
        return 1;
    }

    return 0;
}

//---------------------------------------------------------------------------//
// end of tstfpuTrap.cc
//---------------------------------------------------------------------------//
