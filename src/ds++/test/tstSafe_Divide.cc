//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSafe_Divide.cc
 * \author Mike Buksas
 * \date   Tue Jun 21 16:02:52 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Safe_Divide.hh"
#include "../Assert.hh"
#include "../Release.hh"
#include "ds_test.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test()
{

    double max = numeric_limits<double>::max();

    if (safe_pos_divide (1.0, 0.0) != max) ITFAILS;

    if (safe_pos_divide (10.0, 5.0) != 2.0) ITFAILS;
    
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
	std::cout << "ERROR: While testing tstSafe_Divide, " 
		  << err.what()
		  << std::endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstSafe_Divide, " 
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
	std::cout << "**** tstSafe_Divide Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;

    std::cout << "Done testing tstSafe_Divide "
              << std::endl;
    
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSafe_Divide.cc
//---------------------------------------------------------------------------//
