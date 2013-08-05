//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstProduct_Chebyshev_Legendre.cc
 * \author James S. Warsa
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

#include "../Product_Chebyshev_Legendre.hh"

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
            // SN order = 8, Azimuthal order = 12

            Product_Chebyshev_Legendre quadrature(8, 12); 
            if (quadrature.sn_order()!=8)
            {
                ut.failure("NOT correct SN order");
            }
            quadrature_test(ut, quadrature);
        }
        {
            // SN order = 8, Azimuthal order = 12

            Product_Chebyshev_Legendre quadrature(8, 12, 1, 2);
            quadrature_test(ut, quadrature);
        }

        Product_Chebyshev_Legendre quadrature4(4,8); 
        quadrature_integration_test(ut, quadrature4);
        Product_Chebyshev_Legendre quadrature8(8,16); 
        quadrature_integration_test(ut, quadrature8);
        Product_Chebyshev_Legendre quadrature10(10,20); 
        quadrature_integration_test(ut, quadrature10);
        Product_Chebyshev_Legendre quadrature12(12,24); 
        quadrature_integration_test(ut, quadrature12);
        Product_Chebyshev_Legendre quadrature16(16,32); 
        quadrature_integration_test(ut, quadrature16);

    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstProductGaussLegendre, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstProductGaussLegendre, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstProduct_Chebyshev_Legendre.cc
//---------------------------------------------------------------------------//
