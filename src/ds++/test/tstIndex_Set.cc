//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstIndex_Set.cc
 * \author Mike Buksas
 * \date   Thu Feb  2 13:46:36 2006
 * \brief  
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Assert.hh"
#include "../Release.hh"
#include "ds_test.hh"

#include "../Index_Set.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_index_set()
{

    unsigned dimensions[] = {3,4,5};
    
    Index_Set<3,1> box(dimensions);

    // Check the sizes and ranges of each dimensions
    if (box.get_size()     != 60) ITFAILS;
    if (box.min_of_index() !=  1) ITFAILS;
    if (box.max_of_index() != 60) ITFAILS;

    if (box.get_size(0)     != 3) ITFAILS;
    if (box.min_of_index(0) != 1) ITFAILS;
    if (box.max_of_index(0) != 3) ITFAILS;

    if (box.get_size(1)     != 4) ITFAILS;
    if (box.min_of_index(1) != 1) ITFAILS;
    if (box.max_of_index(1) != 4) ITFAILS;

    if (box.get_size(2)     != 5) ITFAILS;
    if (box.min_of_index(2) != 1) ITFAILS;
    if (box.max_of_index(2) != 5) ITFAILS;


    // Test for indices in the total range and the range of each
    // dimension.
    if (box.index_in_range(0))   ITFAILS;
    if (!box.index_in_range(1))  ITFAILS;
    if (!box.index_in_range(60)) ITFAILS;
    if (box.index_in_range(61))  ITFAILS;

    if (box.index_in_range(0,0))  ITFAILS;
    if (!box.index_in_range(1,0)) ITFAILS;
    if (!box.index_in_range(3,0)) ITFAILS;
    if (box.index_in_range(4,0))  ITFAILS;

    if (box.index_in_range(0,0))  ITFAILS;
    if (!box.index_in_range(1,0)) ITFAILS;
    if (!box.index_in_range(3,0)) ITFAILS;
    if (box.index_in_range(4,0))  ITFAILS;


    // Test the functions for vetting direction and dimension arguments.
    if (box.direction_okay(0))  ITFAILS;
    if (!box.direction_okay(1)) ITFAILS;
    if (!box.direction_okay(6)) ITFAILS;
    if (box.direction_okay(7))  ITFAILS;

    if (box.dimension_okay(-1)) ITFAILS;
    if (!box.dimension_okay(0)) ITFAILS;
    if (!box.dimension_okay(2)) ITFAILS;
    if (box.dimension_okay(3))  ITFAILS;


    // Resize the object and repeat some of the tests:

    // Make a uniform array 10x10x10
    box.set_size(10);

    if (box.get_size() != 1000) ITFAILS;

    dimensions[0] = 10;
    dimensions[1] = 2;
    dimensions[2] = 5;

    box.set_size(dimensions);

    if (box.get_size() != 100) ITFAILS;
        
        
        
    



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
        test_index_set();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstIndex_Set, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstIndex_Set, " 
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
        std::cout << "**** tstIndex_Set Test: PASSED"
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;
    

    std::cout << "Done testing tstIndex_Set"
              << std::endl;
    

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstIndex_Set.cc
//---------------------------------------------------------------------------//
