//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSurface_Descriptor.cc
 * \author Mike Buksas
 * \date   Mon Aug 18 10:23:25 2003
 * \brief  Unit test executable for Surface_Descriptor
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
#include "mc_test.hh"
#include "../Surface_Descriptor.hh"

using namespace std;
using namespace rtt_mc;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void packing_test()
{

    double data[2] = {0.0, 1.0};
    vector<double> surface_data(data, data+2);

    Surface_Descriptor sd(0, surface_data);

    // Test the copy constructor and inequality operator
    Surface_Descriptor sd2(sd);
    if (sd != sd2) ITFAILS;

    // Pack that puppy
    vector<char> packed_sd = sd.pack();

    // Unpack and test for equality
    Surface_Descriptor sd_unpacked(packed_sd);
    if (sd != sd_unpacked) ITFAILS;

}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_mc::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	packing_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstSurface_Descriptor, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_mc_test::passed) 
    {
        std::cout << "**** tstSurface_Descriptor Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstSurface_Descriptor." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSurface_Descriptor.cc
//---------------------------------------------------------------------------//
