//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstSurface_tracker.cc
 * \author Mike Buksas
 * \date   Thu Jun 19 11:54:14 2003
 * \brief  Test execuatble for Surface_tracker class.
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

#include "mc/Sphere.hh"
#include "mc/Surface.hh"
#include "../Surface_tracker.hh"
#include <vector>

using namespace std;
using namespace rtt_mc;
using namespace rtt_imc;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

SP<Surface_tracker> build_test_surface_tracker()
{

    vector<SP<Surface> > surfaces;

    surfaces.push_back( SP<Sphere>( new Sphere( 0.0, 1.0) ) );

    surfaces.push_back( SP<Sphere>( new Sphere( 1.0, 1.0) ) );

    surfaces.push_back( SP<Sphere>( new Sphere(-1.0, 3.0) ) );

    SP<Surface_tracker> tracker(new Surface_tracker(surfaces));

    return tracker;

}


void test_initial_status()
{

    SP<Surface_tracker> tracker ( build_test_surface_tracker() );
    vector<double> position(3), direction(3);

    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;

    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;

    tracker->initialize_status(position, direction);

    if ( tracker->get_inside(0) != true) ITFAILS;
    if ( tracker->get_inside(1) != true) ITFAILS;
    if ( tracker->get_inside(2) != true) ITFAILS;

    direction[2] = -1.0;

    tracker->initialize_status(position, direction);

    if ( tracker->get_inside(0) != true)  ITFAILS;
    if ( tracker->get_inside(1) != false) ITFAILS;
    if ( tracker->get_inside(2) != true)  ITFAILS;
    

}

void test_streaming()
{

    SP<Surface_tracker> tracker ( build_test_surface_tracker() );

    vector<double> position(3);
    vector<double> direction(3);

    position[0] =  0.0;
    position[1] =  0.0;
    position[2] = -1.0;

    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;

    tracker->initialize_status(position, direction);

    const double distance = 5.0;
    const double ew = 1.0;
    const double sigma = 1.0;

    tracker->tally_crossings(position, direction, distance, ew, sigma);

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

	test_initial_status();

	test_streaming();

    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstSurface_tracker, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_imc_test::passed) 
    {
        std::cout << "**** tstSurface_tracker Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstSurface_tracker." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSurface_tracker.cc
//---------------------------------------------------------------------------//
