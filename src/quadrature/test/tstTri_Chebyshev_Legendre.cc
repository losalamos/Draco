//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstTri_Chebyshev_Legendre.cc
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

#include "../Tri_Chebyshev_Legendre.hh"

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
        {
            Tri_Chebyshev_Legendre quadrature(8); // SN order = 8
            if (quadrature.sn_order()!=8)
            {
                ut.failure("NOT correct SN order");
            }
            quadrature_test(ut, quadrature);
        }
        {
            Tri_Chebyshev_Legendre quadrature(8, 1, 2); // SN order = 8, mu=1, eta=2
            quadrature_test(ut, quadrature);
        }
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstGaussLegendre, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstGaussLegendre, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstTri_Chebyshev_Legendre.cc
//---------------------------------------------------------------------------//
