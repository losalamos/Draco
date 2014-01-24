//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/test/tstPower.cc
 * \author Mike Buksas
 * \date   Mon Jul 24 13:47:58 2006
 * \brief  
 * \note   Copyright © 2006-2014 Los Alamos National Security, LLC
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
#include "ds++/Release.hh"

#include "../Power.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test( rtt_dsxx::UnitTest &ut )
{
    
    using rtt_sf::Power;

    // Power of an integer:

    if (Power<12>(2) != 4096)
        ut.failure("2^12 == 4096 in integers failed");
    else
        ut.passes("");

    // Floating point bases:
    
    if (!soft_equiv(Power<4> (2.0), 16.0) )
        ut.failure("2.0^4 = 16.0 in float failed.");
    else
        ut.passes("");

    if (!soft_equiv(Power<17>(1.1), 5.054470284992945) )
        ut.failure("1.1^17 failed.");
    else
        ut.passes("");

    if (Power<0>(1) != 1)
        ut.failure("1^0 in int failed.");
    else
        ut.passes("");

    if (!soft_equiv(Power<0>(1.0), 1.0))
        ut.failure("1.0^0 in float failed.");
    else
        ut.passes("");

}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " 
                 << rtt_dsxx::release() 
                 << endl;
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );

        test(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstPower, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstPower, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        return 1;
    }


    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstPower.cc
//---------------------------------------------------------------------------//
