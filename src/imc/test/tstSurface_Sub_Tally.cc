//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstSurface_Sub_Tally.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 15:52:59 2003
 * \brief  Unit test for Surface_Sub_Tally
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "imc_test.hh" 

#include "ds++/Soft_Equivalence.hh"
#include "../Surface_Sub_Tally.hh"
#include "../Azimuthal_Mesh.hh" 
#include "ds++/SP.hh"

using namespace rtt_imc;
using namespace rtt_dsxx;
using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test()
{

    double c[3] = {-0.5, 0.0, 0.5};
    vector<double> cosines(c, c+3);

    SP<Azimuthal_Mesh> az_mesh ( new Azimuthal_Mesh(cosines) );

    Surface_Sub_Tally surface_tally(az_mesh, 1);

    if (surface_tally.get_number_surfaces() != 1) ITFAILS;
    if (surface_tally.get_outward_weight_tally(1).size() != 4) ITFAILS;
    if (surface_tally.get_mesh_size() != 4 ) ITFAILS;

    vector<double> direction(3);
    direction[0] = std::sqrt(2.0)/2.0;
    direction[1] = 0.0;
    direction[2] = std::sqrt(2.0)/2.0;

    surface_tally.add_to_tally(1, direction, true, 1.0);

    if ( !soft_equiv( surface_tally.get_outward_weight_tally(1)[3], 1.0 ) ) ITFAILS;

    direction[0] = std::sqrt( 1.0 - 0.1 * 0.1);
    direction[1] =  0.0;
    direction[2] = -0.1;

    surface_tally.add_to_tally(1, direction, false, 0.5);

    if (!soft_equiv( surface_tally.weight(1,false,2), 0.5 ) ) ITFAILS;
    if (surface_tally.crossings(1,false,2) != 1) ITFAILS;

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
	std::cout << "While testing tstSurface_Sub_Tally, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_imc_test::passed) 
    {
        std::cout << "**** tstSurface_Sub_Tally Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstSurface_Sub_Tally." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSurface_Sub_Tally.cc
//---------------------------------------------------------------------------//
