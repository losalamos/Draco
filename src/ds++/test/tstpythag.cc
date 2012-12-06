//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   dsxx/test/tstpythag.cc
 * \author Kent Budge
 * \date   Mon Aug  9 14:45:55 2004
 * \brief  Test the pythag function
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

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
    {
	ut.passes("pythag correct");
    }
    else
    {
	ut.failure("pythag NOT correct");
    }
    if (soft_equiv(pythag(4.0e307, 3.0e307), 5.0e307))
    {
	ut.passes("pythag correct");
    }
    else
    {
	ut.failure("pythag NOT correct");
    }
    if (soft_equiv(pythag(0.0, 0.0), 0.0))
    {
	ut.passes("pythag correct");
    }
    else
    {
	ut.failure("pythag NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
        tstpythag( ut );
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown." << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstpythag.cc
//---------------------------------------------------------------------------//
