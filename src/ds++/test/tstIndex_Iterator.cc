//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstIndex_Iterator.cc
 * \author Mike Buksas
 * \date   Wed Feb  1 08:58:48 2006
 * \brief  Unit test for Index_Iterator
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
#include "../Index_Iterator.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_index_iterator()
{
    
    unsigned dimensions[] = {3,4,5};
    

    Index_Converter<3,1> box(dimensions);

    Index_Iterator<3,1> it(box);

    if (it.get_index()  != 1) ITFAILS;
    if (it.get_index(0) != 1) ITFAILS;
    if (it.get_index(1) != 1) ITFAILS;
    if (it.get_index(2) != 1) ITFAILS;
    if (!it.is_in_range()) ITFAILS;

    ++it;

    if (it.get_index()  != 2) ITFAILS;
    if (it.get_index(0) != 2) ITFAILS;
    if (it.get_index(1) != 1) ITFAILS;
    if (it.get_index(2) != 1) ITFAILS;
    if (!it.is_in_range()) ITFAILS;

    --it; --it;

    if (it.is_in_range()) ITFAILS;

}

void test_looping()
{

    unsigned dimensions[] = {3,4,5};

    Index_Converter<3,1> box(dimensions);

    int index = 1;
    for (Index_Iterator<3,1> it(box); it.is_in_range(); ++it)
    {
        const int it_index = it.get_index();
        
        // Check the returned index against a manual count.
        if (it_index != index++) ITFAILS;
        
        // Check the first and last index directly.
        if ( (it_index-1) %  3 + 1 != it.get_index(0) ) ITFAILS;
        if ( (it_index-1) / 12 + 1 != it.get_index(2) ) ITFAILS;

        
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
        test_index_iterator();

        test_looping();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstIndex_Iterator, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstIndex_Iterator, " 
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
        std::cout << "**** tstIndex_Iterator Test: PASSED"
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;
    

    std::cout << "Done testing tstIndex_Iterator"
              << std::endl;
    

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstIndex_Iterator.cc
//---------------------------------------------------------------------------//
