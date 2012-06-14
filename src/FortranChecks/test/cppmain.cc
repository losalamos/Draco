//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FortranCheck/test/cppmain.cc
 * \author Kelly Thompson
 * \date   Tuesday, Jun 12, 2012, 16:03 pm
 * \brief  Test C++ main linking a Fortran library.
 * \note   Copyright (c) 2012 Los Alamos National Security, LLC
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <iostream>

using namespace std;

// forward declaration of f90 functions
extern "C" void sub1(double alpha, size_t *numPass, size_t *numFail);

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        // >>> UNIT TESTS
        double alpha = 1.0;
        size_t np(ut.numPasses);
        size_t nf(ut.numFails);
        // Call fortran subroutine
        sub1(alpha,&np,&nf);
        ut.numPasses = np;
        ut.numFails = nf;            
        std::cout << ut.numPasses << " " << ut.numFails << std::endl;
    }
    catch (rtt_dsxx::assertion &assert)
    {
        cout << "While testing cppmain, " << assert.what() << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing cppmain, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
// end of cppmain.cc
//---------------------------------------------------------------------------//
