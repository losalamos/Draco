//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/tstF4.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:57:47 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "ds++/Release.hh"
#include "../F4.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_sf;
using rtt_units::PI;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstF4( UnitTest & ut )
{
    double f3 = F4(-10.0);   
    if (soft_equiv(f3, 4*3*2*exp(-10.0)*(1-exp(-10.0)/32.0), 2e-6))
    {
	ut.passes("correct F4 for -20.0");
    }
    else
    {
	ut.failure("NOT correct F4 for -20.0");
    }
    f3 = F4(1000.0);   
    if (soft_equiv(f3,
                   pow(1000.0, 5.0)/5.0
                   + PI*PI*4*pow(1000.0, 3.0)/6.0, 1.0e-10))
    {
	ut.passes("correct F4 for 1000.0");
    }
    else
    {
	ut.failure("NOT correct F4 for 1000.0");
    }

    ofstream out("debug.dat");
    for (double eta = -10; eta<20; eta += 0.1)
    {
        out << eta << ' ' << F4(eta) << endl;
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstF4(ut);
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
//                        end of tstF4.cc
//---------------------------------------------------------------------------//
