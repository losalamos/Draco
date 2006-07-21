//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstPower.cc
 * \author Mike Buksas
 * \date   Thu Jul 20 17:31:36 2006
 * \brief  
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Assert.hh"
#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "ds_test.hh"

#include "../Power.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test()
{

    if (Power<4>(2.0) != 16.0) ITFAILS;
    if (Power<4>(10.0) != 10000.0) ITFAILS;

    if (Power<6> (2.0)  != 64.0) ITFAILS;
    if (Power<12>(2.0) != 4096.0) ITFAILS;

    if (!soft_equiv(Power<17>(2.5), 5820766.0913467407)) ITFAILS;


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
        test();

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
		  << "An unknown exception was thrown"
                  << std::endl;
        return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
              << std::endl;
    if (rtt_ds_test::passed) 
    {
        std::cout << "**** tstPower Test: PASSED "
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;

    std::cout << "Done testing tstPower "
              << std::endl;
    

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstPower.cc
//---------------------------------------------------------------------------//
