//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstabs.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  
 * \note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../abs.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstabs(UnitTest &ut)
{
    using rtt_dsxx::abs;
    
    if (abs(-2.2)==abs(2.2))
    {
        ut.passes("Correctly calculated abs(double)");
    }
    else
    {
        ut.failure("Did NOT correctly calculate abs(double)");
    }
    if (abs(-2.2f)==abs(2.2f))
    {
        ut.passes("Correctly calculated abs(float)");
    }
    else
    {
        ut.failure("Did NOT correctly calculate abs(float)");
    }
    if (abs(-2)==abs(2))
    {
        ut.passes("Correctly calculated abs(int)");
    }
    else
    {
        ut.failure("Did NOT correctly calculate abs(int)");
    }
    if (abs(-2L)==abs(2L))
    {
        ut.passes("Correctly calculated abs(long)");
    }
    else
    {
        ut.failure("Did NOT correctly calculate abs(long)");
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstabs(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstabs.cc
//---------------------------------------------------------------------------//
