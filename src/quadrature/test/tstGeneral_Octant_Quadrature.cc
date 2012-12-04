//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstGeneral_Octant_Quadrature.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id: template_test.cc 5830 2011-05-05 19:43:43Z kellyt $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"

#include "quadrature_test.hh"

#include "../General_Octant_Quadrature.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        double const V = 1/sqrt(3.0);
        vector<double> mu(1, V), eta(1, V), xi(1, V), wt(1, 1.0);
        General_Octant_Quadrature quadrature(mu,
                                             eta,
                                             xi,
                                             wt,
                                             2,
                                             Quadrature::TRIANGLE_QUADRATURE);

        quadrature_test(ut, quadrature);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstGereral_Octant_Quadrature, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstGereral_Octant_Quadrature, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstGeneral_Octant_Quadrature.cc
//---------------------------------------------------------------------------//
