//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/test/tstqr_unpack.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "../Release.hh"
#include "../qr_unpack.hh"
#include "../qrdcmp.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstqr_unpack(UnitTest &ut)
{
    vector<double> A(2*2);
    A[0+2*0] = 2.;
    A[0+2*1] = 3.;
    A[1+2*0] = 1.;
    A[1+2*1] = 5.;

    vector<double> C, D;

    qrdcmp(A, 2, C, D);

    vector<double> Qt;
    vector<double> R = A;    
    qr_unpack(R, 2, C, D, Qt);

    if (soft_equiv(R[0+2*0], D[0]))
    {
	ut.passes("R[0+2*0] is correct");
    }
    else
    {
	ut.failure("R[0+2*0] is NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstqr_unpack(ut);
    }
    catch (exception &err)
    {
	cout << "ERROR: While testing tstqr_unpack, " << err.what() << endl;
	return 1;
    }
    catch( ... )
    {
	cout << "ERROR: While testing tstqr_unpack, " 
             << "An unknown exception was thrown." << endl;
	return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstqr_unpack.cc
//---------------------------------------------------------------------------//
