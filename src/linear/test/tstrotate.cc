//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/test/tstrotate.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright (C) 2004-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include "../rotate.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstrotate(UnitTest &ut)
{
    vector<double> A(3*3, 0.0);
    for (unsigned i=0; i<3; i++)
    {
	A[i+3*i] = 1.0;
    }
    vector<double> B = A;

    rotate(A, B, 3, 1, 0.5, 0.5);

    if (soft_equiv(A[0+3*0], 1.0))
    {
	ut.passes("A[0+3*0] is correct");
    }
    else
    {
	ut.failure("A[0+3*0] is NOT correct");
    }
    if (soft_equiv(A[1+3*1], sqrt(0.5)))
    {
	ut.passes("A[1+3*1] is correct");
    }
    else
    {
	ut.failure("A[1+3*1] is NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstrotate(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstrotate, " << err.what() << endl;
	return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstrotate, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstrotate.cc
//---------------------------------------------------------------------------//
