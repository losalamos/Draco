//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/tstF_eta_inv.cc
 * \author Kent Budge
 * \date   Mon Sep 20 14:55:09 2004
 * \brief  Test the F_eta_inv function
 * \note   Copyright (C) 2004-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <cmath>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "../F_eta.hh"
#include "../F_eta_inv.hh"

using namespace std;
using namespace rtt_sf;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstF_eta_inv( UnitTest & ut )
{
    const double C_KB = 1.38066e-16;  // Boltzmann's constant in cgs
    const double C_ME = 9.108e-28;    // Electron mass in grams
    const double C_C = 2.99792e10;    // Speed of light in cm/sec

    const unsigned ntests = 19;
    double eta[ntests] = {
        -70.0, -5.0, -0.694, -0.693, 0.0, 5.0, 10.0, 20.0, 50.0, 100.0,
        800.0,
        -50.0, -5.0, 0.0, 5.0, 
        -20.0, -5.0, 5.0, 50.0
    };
    double T[ntests] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1.0e9, 1.0e9, 1.0e9, 1.0e9,
        500.0e9, 500.0e9, 500.0e9, 500.0e9
    };

    for (unsigned i=0; i<ntests; i++)
    {
        double gamma = T[i]*(C_KB/(C_ME*C_C*C_C));
        double reta = gamma*eta[i];
        double rreta = F_eta_inv(F_eta(reta, gamma), gamma);
        if (fabs((rreta-reta)/gamma)>1.0e-5)
        {
            ut.failure("F_eta_inv FAILED");
        }
        else
        {
            ut.passes("F_eta_inv passed");
        }
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try { tstF_eta_inv(ut); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstF_eta_inv.cc
//---------------------------------------------------------------------------//
