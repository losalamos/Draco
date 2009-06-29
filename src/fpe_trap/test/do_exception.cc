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

#include <fstream>
#include <cmath>
#include <ds++/Assert.hh>
#include "../fpe_trap.hh"

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
//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    Insist(argc == 2, "Wrong number of args.");

    std::ofstream f;

    f.open("output.dat");

    if ( rtt_fpe_trap::enable_fpe() )
    {
        // Platform supported.
        f << "supported" << endl;
    }
    else
    {
        // Platform not supported.
        f << "unsupported\n";
        f.close();

        return 0;
    }

    int test;

    sscanf(argv[1], "%d", &test);

    Insist(test >= 0 && test <= 3, "Bad test value.");

    double zero = 0.0; // for double division by zero
    double neg = -1.0; // for sqrt(-1.0)
    double result;

    // Certain tests may be optimized away by the compiler, by
    // recogonizing the constants set above and precomputing the
    // results below.  So do something here to hopefully avoid this.
    // This tricks the optimizer, at least for gnu and KCC.

    if ( test < -100 )
    { // this should never happen
	Insist(0, "Something is very wrong.");
	zero = neg = 1.0; // trick the optimizer?
    }

    switch ( test ) {
    case 0:
	f << "should_work" << endl;
	result = 1.0 + zero + sqrt(-neg);
	f << "result = " << result << endl;
	break;
    case 1:
	f << "div_by_zero" << endl;
	result = 1.0 / zero; // should fail here
	f << "result = " << result << endl;
	break;
    case 2:
	f << "sqrt(-1.0)" << endl;
	result = sqrt(neg); // should fail here
	f << "result = " << result << endl;
	break;
    case 3: {
	f << "overflow" << endl;
	result = 2.0;
	for ( int i = 0; i < 100; i++ ) {
	    result = exp(result); // should fail at some i
	}
	f << "result = " << result << endl;
	break;
    }
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
// end of tstfpuTrap.cc
//---------------------------------------------------------------------------//
