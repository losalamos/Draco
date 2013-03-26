//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSQR.cc
 * \author Kent Budge
 * \date   Tue Jul  6 08:54:38 2004
 * \brief  Test the SQR function.
 * \note   Copyright (C)  2004-2013 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../square.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try
    {
	if (soft_equiv(square(3.0), 9.0))
	    ut.passes("square function returned correct double");
	else
	    ut.failure("square function did NOT return correct double.");
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstSQR.cc
//---------------------------------------------------------------------------//
