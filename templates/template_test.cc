//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "<spkg>.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_<pkg>::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing <class>, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_<spkg>::passed) 
    {
        std::cout << "**** <class> Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing <class>." << std::endl;
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
