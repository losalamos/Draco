//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstcube.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
 * \brief  
 * \note   Copyright (C) 2004-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../cube.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try
    {
	if (soft_equiv(cube(2.0), 8.0))
	{
	    ut.passes("square function returned correct double");
	}
	else
	{
	    ut.failure("square function did NOT return correct double.");
	}
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of testcube.cc
//---------------------------------------------------------------------------//
