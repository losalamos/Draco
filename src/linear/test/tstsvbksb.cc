//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/test/tstsvbksb.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright (C) 2004-2014 Los Alamos National Security, LLC.
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
#include "../svdcmp.hh"
#include "../svbksb.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstsvbksb(UnitTest &ut)
{
    vector<double> U(4), W(2), V(4);
    U[0+2*0] = 2.;
    U[0+2*1] = 3.;
    U[1+2*0] = 1.;
    U[1+2*1] = 5.;

    svdcmp(U, 2, 2, W, V);

    vector<double> b(2), x;
    b[0] = 1.0;
    b[1] = 2.0;

    svbksb(U, W, V, 2, 2, b, x);

    if (soft_equiv(b[0], x[0]*2 + x[1]*3))
    {
	ut.passes("0 is correct");
    }
    else
    {
	ut.failure("0 is NOT correct");
    }
    if (soft_equiv(b[1], x[0]*1 + x[1]*5))
    {
	ut.passes("0 is correct");
    }
    else
    {
	ut.failure("0 is NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstsvbksb(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstsvbksb, " << err.what() << endl;
	return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstsvbksb, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstsvbksb.cc
//---------------------------------------------------------------------------//
