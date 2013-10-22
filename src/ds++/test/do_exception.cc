//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/test/do_exception.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:33:59 2005
 * \brief  Does a floating-point exception.
 * \note   Copyright (C) 2005-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../fpe_trap.hh"
#include "../Assert.hh"
#include <fstream>
#include <cmath>
#include <iostream>

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

void run_test(int /*argc*/, char *argv[])
{
    std::ofstream fout;
    fout.open("output.dat");

    bool const abortWithInsist(true);
    if ( rtt_dsxx::fpe_trap(abortWithInsist).enable() )
    {
        // Platform supported.
        fout << "- fpe_trap: This platform is supported" << endl;
        cout << "- fpe_trap: This platform is supported" << endl;
    }
    else
    {
        // Platform not supported.
        fout << "- fpe_trap: This platform is not supported\n";
        fout.close();
        return;
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
            fout << "- Case zero: this operation should not throw a SIGFPE."
                 << " The result should be 2." << endl;
            result = 1.0 + zero + sqrt(-neg);
            fout << "  result = " << result << endl;
            break;
        case 1:
            fout << "- Trying a div_by_zero operation" << endl;
            cout << "- Trying a div_by_zero operation" << endl;
            result = 1.0 / zero; // should fail here
            cout << "  result = " << 1.0*result << endl;
            fout << "  result = " << 1.0*result << endl;
            break;
        case 2:
            fout << "- Trying to evaluate sqrt(-1.0)" << endl;
            result = sqrt(neg); // should fail here
            fout << "  result = " << result << endl;
            break;
        case 3:
        {
            fout << "- Trying to cause an overflow condition" << endl;
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
            fout << "  result = " << result << endl;
            break;
        }
    }
    // close the log file.
    fout.close();
    return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    Insist(argc == 2, "Wrong number of args.");
    try
    {
        run_test(argc,argv);
    }
    catch (exception &err)
    {
        if ( rtt_dsxx::fpe_trap().enable() )
        {
            // keyword 'signal' shows up as a failure when processed by
            // test_filter.py.  To avoid this, we do not print the err.what()
            // message. 
            cout << "While running " << argv[0] << ", " 
                 << "a SIGFPE was successfully caught.\n\t"
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
// end of do_exception.cc
//---------------------------------------------------------------------------//
