//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/test/tstsvdcmp.cc
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
#include "../svdcmp.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstsvdcmp(UnitTest &ut)
{
    vector<double> U(6), W, V;
    U[0+3*0] = 2.;
    U[0+3*1] = 3.;
    U[1+3*0] = 1.;
    U[1+3*1] = 5.;
    U[2+3*0] = 4.;
    U[2+3*1] = 4;

    svdcmp(U, 3, 2, W, V);

    // Compute U*W*Tr(V) to verify

    double WV[9];
    WV[0+2*0] = W[0]*V[0+2*0];
    WV[0+2*1] = W[0]*V[1+2*0];
    WV[1+2*0] = W[1]*V[0+2*1];
    WV[1+2*1] = W[1]*V[1+2*1];

    double UWV[6];
    UWV[0+3*0] = U[0+3*0]*WV[0+2*0] + U[0+3*1]*WV[1+2*0];
    UWV[0+3*1] = U[0+3*0]*WV[0+2*1] + U[0+3*1]*WV[1+2*1];
    UWV[1+3*0] = U[1+3*0]*WV[0+2*0] + U[1+3*1]*WV[1+2*0];
    UWV[1+3*1] = U[1+3*0]*WV[0+2*1] + U[1+3*1]*WV[1+2*1];
    UWV[2+3*0] = U[2+3*0]*WV[0+2*0] + U[2+3*1]*WV[1+2*0];
    UWV[2+3*1] = U[2+3*0]*WV[0+2*1] + U[2+3*1]*WV[1+2*1];

    if (soft_equiv(UWV[0+3*0], 2.0))
    {
	ut.passes("0,0 is correct");
    }
    else
    {
	ut.failure("0,0 is NOT correct");
    }
    if (soft_equiv(UWV[0+3*1], 3.0))
    {
	ut.passes("0,1 is correct");
    }
    else
    {
	ut.failure("0,1 is NOT correct");
    }
    if (soft_equiv(UWV[1+3*0], 1.0))
    {
	ut.passes("1,0 is correct");
    }
    else
    {
	ut.failure("1,0 is NOT correct");
    }
    if (soft_equiv(UWV[1+3*1], 5.0))
    {
	ut.passes("1,1 is correct");
    }
    else
    {
	ut.failure("1,1 is NOT correct");
    }
    if (soft_equiv(UWV[2+3*1], 4.0))
    {
	ut.passes("2,1 is correct");
    }
    else
    {
	ut.failure("2,1 is NOT correct");
    }
    if (soft_equiv(UWV[2+3*1], 4.0))
    {
	ut.passes("2,1 is correct");
    }
    else
    {
	ut.failure("2,1 is NOT correct");
    }

    // Now decompose transpose, to exercise different code paths

    U[0+2*0] = 2.;
    U[1+2*0] = 3.;
    U[0+2*1] = 1.;
    U[1+2*1] = 5.;
    U[0+2*2] = 4.;
    U[1+2*2] = 4;

    svdcmp(U, 2, 3, W, V);

    // Compute U*W*Tr(V) to verify

    WV[0+3*0] = W[0]*V[0+3*0];
    WV[0+3*1] = W[0]*V[1+3*0];
    WV[0+3*2] = W[0]*V[2+3*0];
    WV[1+3*0] = W[1]*V[0+3*1];
    WV[1+3*1] = W[1]*V[1+3*1];
    WV[1+3*2] = W[1]*V[2+3*1];
    WV[2+3*0] = W[2]*V[0+3*2];
    WV[2+3*1] = W[2]*V[1+3*2];
    WV[2+3*2] = W[2]*V[2+3*2];

    UWV[0+2*0] = U[0+2*0]*WV[0+3*0] + U[0+2*1]*WV[1+3*0] + U[0+2*2]*WV[2+3*0];
    UWV[0+2*1] = U[0+2*0]*WV[0+3*1] + U[0+2*1]*WV[1+3*1] + U[0+2*2]*WV[2+3*1];
    UWV[0+2*2] = U[0+2*0]*WV[0+3*2] + U[0+2*1]*WV[1+3*2] + U[0+2*2]*WV[2+3*2];
    UWV[1+2*0] = U[1+2*0]*WV[0+3*0] + U[1+2*1]*WV[1+3*0] + U[1+2*2]*WV[2+3*0];
    UWV[1+2*1] = U[1+2*0]*WV[0+3*1] + U[1+2*1]*WV[1+3*1] + U[1+2*2]*WV[2+3*1];
    UWV[1+2*2] = U[1+2*0]*WV[0+3*2] + U[1+2*1]*WV[1+3*2] + U[1+2*2]*WV[2+3*2];

    if (soft_equiv(UWV[0+2*0], 2.0))
    {
	ut.passes("0,0 is correct");
    }
    else
    {
	ut.failure("0,0 is NOT correct");
    }
    if (soft_equiv(UWV[0+2*1], 1.0))
    {
	ut.passes("0,1 is correct");
    }
    else
    {
	ut.failure("0,1 is NOT correct");
    }
    if (soft_equiv(UWV[0+2*2], 4.0))
    {
	ut.passes("0,2 is correct");
    }
    else
    {
	ut.failure("0,2 is NOT correct");
    }
    if (soft_equiv(UWV[1+2*0], 3.0))
    {
	ut.passes("1,0 is correct");
    }
    else
    {
	ut.failure("1,0 is NOT correct");
    }
    if (soft_equiv(UWV[1+2*1], 5.0))
    {
	ut.passes("1,1 is correct");
    }
    else
    {
	ut.failure("1,1 is NOT correct");
    }
    if (soft_equiv(UWV[1+2*2], 4.0))
    {
	ut.passes("1,2 is correct");
    }
    else
    {
	ut.failure("1,2 is NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstsvdcmp(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstsvdcmp, " << err.what() << endl;
	return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstsvdcmp, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstsvdcmp.cc
//---------------------------------------------------------------------------//
