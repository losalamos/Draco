//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstAzimuthal_Mesh.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 13:59:30 2003
 * \brief  Unit test for the aziumthal mesh class
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "../Release.hh"
#include "imc_test.hh" 

#include "../Azimuthal_Mesh.hh"

using namespace rtt_imc;
using namespace rtt_dsxx;
using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test()
{

    double c[6] = {-1.0, -0.5, 0.0, 0.25, 0.5, 1.0};

    vector<double> cosines(c, c+6);

    Azimuthal_Mesh az_mesh(cosines);

    if (az_mesh.size() != 5) ITFAILS;

    if (!soft_equiv(az_mesh.get_lower_cosine(1), -1.0) ) ITFAILS;
    if (!soft_equiv(az_mesh.get_upper_cosine(5),  1.0) ) ITFAILS;

    if (!soft_equiv(az_mesh.get_upper_cosine(2),
		    az_mesh.get_lower_cosine(3) ) ) ITFAILS;


    vector<double> direction(3);

    direction[0] = 0.0; direction[1] = 0.0; direction[2] = 1.0;
    if (!az_mesh.is_in_bin(direction, 5) ) ITFAILS;
    if (az_mesh.find_bin(direction) != 5 ) ITFAILS;

    direction[2] = 0.5;
    if (!az_mesh.is_in_bin(direction, 4) ) ITFAILS;
    if (!az_mesh.is_in_bin(direction, 5) ) ITFAILS;
    if (az_mesh.find_bin(direction) != 4 ) ITFAILS;

    direction[2] = 0.1;
    if (!az_mesh.is_in_bin(direction, 3) ) ITFAILS;
    if (az_mesh.find_bin(direction) != 3 ) ITFAILS;


#if DBC & 2
    try {
	cosines[0] = -0.75;
	Azimuthal_Mesh az_mesh(cosines);
	ITFAILS;
    }
    catch (rtt_dsxx::assertion &ass) 
    {
	std::cerr << "Caught a deliberate assertion." << std::endl;
    }
#endif

#if DBC & 2
    try {
	cosines[2] = 0.0;
	Azimuthal_Mesh az_mesh(cosines);
	ITFAILS;
    }
    catch (rtt_dsxx::assertion &ass) 
    {
	std::cerr << "Caught a deliberate assertion." << std::endl;
    }
#endif

    
#if DBC & 2
    try {
	cosines[3] = 1.0;
	Azimuthal_Mesh az_mesh(cosines);
	ITFAILS;
    }
    catch (rtt_dsxx::assertion &ass) 
    {
	std::cerr << "Caught a deliberate assertion." << std::endl;
    }
#endif




}
    


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_imc::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS 
	test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstAzimuthal_Mesh, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_imc_test::passed) 
    {
        std::cout << "**** tstAzimuthal_Mesh Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstAzimuthal_Mesh." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstAzimuthal_Mesh.cc
//---------------------------------------------------------------------------//
