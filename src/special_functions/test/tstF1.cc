//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/tstF1.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \brief  Unit test for F1 function.
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>

#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "units/PhysicalConstants.hh"

#include "ds++/Release.hh"
#include "../F1.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstF1( UnitTest & ut )
{
    double f1 = F1(-10.0);   
    if (soft_equiv(f1, 1*exp(-10.0)*(1-exp(-10.0)/4.0), 2e-6))
    {
	ut.passes("correct F1 for -10.0");
    }
    else
    {
	ut.failure("NOT correct F1 for -10.0");
    }
    f1 = F1(1000.0);   
    if (soft_equiv(f1,
                   pow(1000.0, 2.0)/2.0
                   + PI*PI*1*pow(1000.0, 0.0)/6.0, 1.0e-10))
    {
	ut.passes("correct F1 for 1000.0");
    }
    else
    {
	ut.failure("NOT correct F1 for 1000.0");
    }

    ofstream out("debug.dat");
    for (double eta = -10; eta<20; eta += 0.1)
    {
        out << eta << ' ' << F1(eta) << endl;
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstF1(ut);
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
//                        end of tstF1.cc
//---------------------------------------------------------------------------//
