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
#include <c4/global.hh>
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
    rtt_c4::initialize(argc, argv);
    
    Insist(argc == 2, "Wrong number of args.");

    bool host = (rtt_c4::node() == 0);

    std::ofstream f;

    if ( host ) f.open("output.dat");

    if ( rtt_fpe_trap::enable_fpe() )
    {
        // Platform supported.
        if ( host ) f << "supported" << endl;
    }
    else
    {
        // Platform not supported.
        if ( host )
        {
            f << "unsupported\n";
            f.close();
        }
        
        rtt_c4::global_barrier();
        rtt_c4::finalize();
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
	if ( host ) f << "should_work" << endl;
        rtt_c4::global_barrier();
	result = 1.0 + zero + sqrt(-neg);
	if ( host ) f << "result = " << result << endl;
	break;
    case 1:
	if ( host ) f << "div_by_zero" << endl;
        rtt_c4::global_barrier();
	result = 1.0 / zero; // should fail here
	if ( host ) f << "result = " << result << endl;
	break;
    case 2:
	if ( host ) f << "sqrt(-1.0)" << endl;
        rtt_c4::global_barrier();
	result = sqrt(neg); // should fail here
	if ( host ) f << "result = " << result << endl;
	break;
    case 3: {
	if ( host ) f << "overflow" << endl;
        rtt_c4::global_barrier();
	result = 2.0;
	for ( int i = 0; i < 100; i++ ) {
	    result = exp(result); // should fail at some i
	}
	if ( host ) f << "result = " << result << endl;
	break;
    }
    }
    
    rtt_c4::global_barrier();
    rtt_c4::finalize();
    return 0;
}

//---------------------------------------------------------------------------//
// end of tstfpuTrap.cc
//---------------------------------------------------------------------------//
