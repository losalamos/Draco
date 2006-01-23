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

#include "ds++/Assert.hh"
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

    Index_Converter<3,1> box_one_based(dimensions);

    if (box_one_based.get_size()  != 60);
    if (box_one_based(dimensions) != 60) ITFAILS;

    int indices[] = {1,1,1};
    if (box_one_based(indices) != 1) ITFAILS;

    indices[0] = 2; indices[1] = 3; indices[2] = 4;
    int one_index = (2-1) + 3*(3-1) + 12*(4-1) + 1;
    if (box_one_based(indices) != one_index) ITFAILS;


    result = box_one_based(one_index);
    if (!std::equal(result.begin(), result.end(), indices)) ITFAILS;



    Index_Converter<3,0> box_zero_based(dimensions);

    if (box_zero_based.get_size() != 60) ITFAILS;

    indices[0] = 0; indices[1] = 0; indices[2] = 0;
    if (box_zero_based(indices) != 0) ITFAILS;

    indices[0] = dimensions[0] - 1;
    indices[1] = dimensions[1] - 1;
    indices[2] = dimensions[2] - 1;
    if (box_zero_based(indices) != 59) ITFAILS;


    result = box_zero_based(59);
    if (!std::equal(result.begin(), result.end(), indices)) ITFAILS;



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
