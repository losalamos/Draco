//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/test/do_exception.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 14:33:59 2005
 * \brief  Does a floating-point exception.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/fpe_trap.hh"
#include "ds++/Assert.hh"
#include "ds++/StackTrace.hh"
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

//---------------------------------------------------------------------------//
// Helper
void fcout( std::string const & msg, std::ofstream & fout )
{
    fout << msg << endl;
    cout << msg << endl;
    return;
}

//---------------------------------------------------------------------------//
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
void run_test(int /*argc*/, char *argv[])
{
    std::ofstream fout;
    fout.open("output.dat");

    bool const abortWithInsist(true);
    rtt_dsxx::fpe_trap fpet(abortWithInsist);
    if ( fpet.enable() )
    {
        // Platform supported.
        fcout("- fpe_trap: This platform is supported",fout);
        if( ! fpet.active() )
            fcout("- fpe_trap: active flag set to false was not expected.",
                  fout );
    }
    else
    {
        // Platform not supported.
        fout << "- fpe_trap: This platform is not supported\n";
        if( fpet.active() )
            fcout("- fpe_trap: active flag set to true was not expected.",
                  fout );
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

    std::ostringstream msg;
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
            fcout( "- Trying a div_by_zero operation", fout );
            result = 1.0 / zero; // should fail here
            msg << "  result = " << 1.0*result;
            fcout( msg.str(), fout );
            break;
        case 2:
            fcout( "- Trying to evaluate sqrt(-1.0)", fout );
            result = std::sqrt(neg); // should fail here
            msg  << "  result = " << result;
            fcout( msg.str(), fout );
            break;
        case 3:
        {
            fcout( "- Trying to cause an overflow condition", fout );
            result = 2.0;
            std::vector<double> data;
            for ( size_t i = 0; i < 100; i++ )
            {
                // should fail at some i
                result = result * result * result * result * result;
                data.push_back(result); // force optimizer to evaluate the above line.
            }
            msg << "  result = " << result << endl;
            fcout( msg.str(), fout );
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
