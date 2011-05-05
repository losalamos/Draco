//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   min/test/tstpowell.cc
 * \author Kent Budge
 * \date   Tue Jul  6 10:00:38 2004
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
#include "units/PhysicalConstants.hh"

#include "ds++/Release.hh"
#include "../powell.hh"
#include "ds++/cube.hh"

using namespace std;
using namespace rtt_min;
using namespace rtt_dsxx;
using rtt_units::PI;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//


unsigned const NP = 2;

static double xf(double const x, vector<double> const &a)
{
    double const xh = sqrt(x);
    double Result = (sqrt(8.0)+a[0]*xh+a[1]*x+x*xh)*x*xh;
    Result = fabs(Result/cube(sqrt(x*x+2*x)) - 1);
    return Result;
}

static double func(vector<double> const &a)
{
    double Result = 0.0;
    for (double x=1.0e-5; x<1000.0; x *= 1.1)
    {
        Result = max(Result, xf(x, a));
    }
    return Result;
}

void tstpowell( UnitTest & ut )
{
    vector<double> p(NP, 0.0);
    vector<double> xi(NP*NP, 0.0);
    for (unsigned i=0; i<NP; ++i)
    {
        xi[i+NP*i] = 1.0;
    }
    unsigned iter = 10000;
    double fret(0);
    double tolerance(1.0e-5);
    
    powell(p, xi, tolerance, iter, fret, func);

    for (unsigned i=0; i<NP; ++i)
        cout << "a[" << i << "] = " << p[i] << endl;

    cout << "Maximum error: " << fret << endl;

    double tmp[2] = {1.34601,4.19265e-09};
    vector<double> expectedSolution(tmp,tmp+2);
    if( rtt_dsxx::soft_equiv( p.begin(), p.end(),
                              expectedSolution.begin(), expectedSolution.end(),
                              tolerance ) )
    {
        ut.passes("Found expected solution.");
    }
    else
    {
        ut.failure("Did not find expected solution.");
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
        tstpowell(ut);
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
//                        end of testpowell.cc
//---------------------------------------------------------------------------//
