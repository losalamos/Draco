//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstIndex_Converter.cc
 * \author Mike Buksas
 * \date   Fri Jan 20 15:53:51 2006
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

#include "../Index_Converter.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_index_converter()
{

    std::vector<int> result(3);

    unsigned dimensions[] = {3,4,5};

    {
        Index_Converter<3,1> box(dimensions);


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

        if (box.get_index(dimensions) != 60) ITFAILS;

        int indices[] = {1,1,1};
        if (box.get_index(indices) != 1) ITFAILS;

        indices[0] = 2; indices[1] = 3; indices[2] = 4;
        int one_index = (2-1) + 3*(3-1) + 12*(4-1) + 1;
        if (box.get_index(indices) != one_index) ITFAILS;


        result = box.get_indices(one_index);
        if (!std::equal(result.begin(), result.end(), indices)) ITFAILS;
    }

    {
        Index_Converter<3,0> box(dimensions);

        if (box.get_size() != 60) ITFAILS;

        int indices[] = {0, 0, 0};
        if (box.get_index(indices) != 0) ITFAILS;

        indices[0] = dimensions[0] - 1;
        indices[1] = dimensions[1] - 1;
        indices[2] = dimensions[2] - 1;
        if (box.get_index(indices) != 59) ITFAILS;

        box.get_indices(59, result.begin());
        if (!std::equal(result.begin(), result.end(), indices)) ITFAILS;

        // Cell 30 has coordinates (0,3,2):

        if (box.get_next_index(30, 1) != -1) ITFAILS;
        if (box.get_next_index(30, 2) != 31) ITFAILS;

        if (box.get_next_index(30, 3) != 27) ITFAILS;
        if (box.get_next_index(30, 4) != 33) ITFAILS;

        if (box.get_next_index(30, 5) != 18) ITFAILS;
        if (box.get_next_index(30, 6) != 42) ITFAILS;
        
        Index_Converter<3,0> copy(box);
        if (copy != box) ITFAILS;

    }


    {
        Index_Converter<5,1> big_box(10);
        if (big_box.get_size(3) != 10)     ITFAILS;
        if (big_box.get_size()  != 100000) ITFAILS;
    }

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
        test_index_converter();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstIndex_Converter, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstIndex_Converter, " 
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
        std::cout << "**** tstIndex_Converter Test: PASSED"
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;
    

    std::cout << "Done testing tstIndex_Converter"
              << std::endl;
    

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstIndex_Converter.cc
//---------------------------------------------------------------------------//
