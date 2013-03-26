//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   dsxx/test/tstpythag.cc
 * \author Kent Budge
 * \date   Mon Aug  9 14:45:55 2004
 * \brief  Test the pythag function
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Soft_Equivalence.hh"
#include "../Release.hh"
#include "../pythag.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstpythag( UnitTest & ut )
{
    if (soft_equiv(pythag(3.0e307, 4.0e307), 5.0e307))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
    if (soft_equiv(pythag(4.0e307, 3.0e307), 5.0e307))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
    if (soft_equiv(pythag(0.0, 0.0), 0.0))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try
    {
        tstpythag( ut );
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstpythag.cc
//---------------------------------------------------------------------------//
