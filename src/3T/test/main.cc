//----------------------------------*-C++-*----------------------------------//
// main.cc
// Geoffrey M. Furnish
// Wed Nov 19 16:59:29 1997
//---------------------------------------------------------------------------//
// @> Main program for running test cases for the 3T solver.
//---------------------------------------------------------------------------//

#include "Test_Prob.hh"

#include "c4/global.hh"

#include <iostream.h>

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    {
    // Introduce a new scope in case the test problem object needs to do MPI
    // work in its dtor.

	SP<Test_Prob> prob = Test_Prob_allocator( argc, argv );
	prob->run();
    }

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
