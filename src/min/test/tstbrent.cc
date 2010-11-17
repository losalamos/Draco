//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   min/test/tstbrent.cc
 * \author Kent G. Budge
 * \date   Tue Nov 16 17:26:03 2010
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "../Release.hh"
#include "../brent.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_min;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

double f(double x)
{
    return cos(x);
}

void tstbrent(UnitTest &ut)
{
    double xmin;
    double root = brent(0.0, 6.28, 2.0, f, 1.0e-12, xmin);

    cout << xmin << " " << root << endl;

    if (soft_equiv(xmin, M_PI))
    {
        ut.passes("correctly found first minimum of cos");
    }
    else
    {
        ut.failure("did NOT correctly find first minimum of cos");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstbrent(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstbrent, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstbrent, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstbrent.cc
//---------------------------------------------------------------------------//
