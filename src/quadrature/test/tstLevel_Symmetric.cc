//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tstLevel_Symmetric.cc
 * \author Kent G. Budge
 * \date   Tue Nov  6 13:08:49 2012
 * \brief  
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC
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

#include "../Level_Symmetric.hh"

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
        for (unsigned N=2; N<18; N+=2)
            // Pathological in curvilinear geometry at N>=18
        {
            cout << "Order " << N << ':' << endl;
            Level_Symmetric quadrature(N);
            if (quadrature.sn_order()!=N)
            {
                ut.failure("NOT correct SN order");
            }
            quadrature_test(ut, quadrature);
        }
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstLevel_Symmetric, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstLevel_Symmetric, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstLevel_Symmetric.cc
//---------------------------------------------------------------------------//
