//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "<spkg>.hh"

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
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing <class>, " 
		  << err.what()
		  << std::endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing <class>, " 
		  << "An unknown exception was thrown."
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_<spkg>::passed) 
    {
        std::cout << "**** <class> Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing <class>." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
