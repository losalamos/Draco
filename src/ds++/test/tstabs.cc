//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstabs.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Assert.hh"
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
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstabs(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstabs, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstabs, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstabs.cc
//---------------------------------------------------------------------------//
