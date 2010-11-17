//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstHomogeneous_New.cc
 * \author Kent Budge
 * \date   Tue Nov 28 09:17:23 2006
 * \brief  test the Homogeneous_New allocator class.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../ScalarUnitTest.hh"
#include "../Assert.hh"
#include "../Release.hh"
#include "../Homogeneous_New.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstHomogeneous_New(UnitTest &ut)
{
    {
        Homogeneous_New allocator(34);
    }
    ut.passes("Construction/destruction did not throw an exception");

    {
        Homogeneous_New allocator(56);
        void *first_object = allocator.allocate();
        cout << "First address: " << first_object << endl;
        void *second_object = allocator.allocate();
        cout << "Second address: " << second_object << endl;
        allocator.deallocate(first_object);
        first_object = allocator.allocate();
        cout << "Reallocated first address: " << first_object << endl;
        void *third_object = allocator.allocate();
        cout << "Third address: " << third_object << endl;
    }

    {
        Homogeneous_New allocator(29);
        allocator.reserve(7);
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            std::cout << argv[0] << ": version " 
                      << rtt_dsxx::release() 
                      << std::endl;
            return 0;
        }

    try
    {
        ScalarUnitTest ut(argc, argv, release);
        tstHomogeneous_New(ut);
        ut.status();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstClass_Parser, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstClass_Parser, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstHomogeneous_New.cc
//---------------------------------------------------------------------------//
