//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstSurface_tracker.cc
 * \author Mike Buksas
 * \date   Thu Jun 19 11:54:14 2003
 * \brief  Test execuatble for Surface_tracker class.
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
#include "imc_test.hh"

#include "ds++/Soft_Equivalence.hh"
#include "mc/Sphere.hh"
#include "mc/Surface.hh"
#include "../Azimuthal_Mesh.hh"
#include "../Surface_tracker.hh"
#include "../Surface_Sub_Tally.hh"
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

Surface_Sub_Tally make_surface_tally()
{

    double c[3] = {-0.5, 0.0, 0.5};
    vector<double> cosines(c, c+3);

    SP<Azimuthal_Mesh> mesh ( new Azimuthal_Mesh(cosines ) );

    return Surface_Sub_Tally(mesh, 3);

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

    double distance = 5.0;
    const double ew = 1.0;
    const double sigma = 1.0;  

    Surface_Sub_Tally tally( make_surface_tally() );
    
    tracker->tally_crossings_implicit_abs(
	position, direction, distance, ew, sigma, tally);

    if (!soft_equiv(tally.get_outward_tally(1)[3], exp(-2.0) ) ) ITFAILS;
    if (!soft_equiv(tally.get_inward_tally (2)[3], exp(-1.0) ) ) ITFAILS;
    if (!soft_equiv(tally.get_outward_tally(2)[3], exp(-3.0) ) ) ITFAILS;
    if (!soft_equiv(tally.get_outward_tally(3)[3], exp(-3.0) ) ) ITFAILS;

    position[0] = -1.0;
    position[1] =  0.0;
    position[2] =  1.0;

    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;

    tracker->initialize_status(position, direction);
    tracker->tally_crossings_implicit_abs(
	position, direction, distance, ew, sigma, tally);

    // Stream
    for (int i=0; i<3; ++i) position[i] += direction[i] * distance;

    // Scatter toward (0,0,-1)
    direction[0] = -2.0/sqrt(5.0);
    direction[1] =  0.0;
    direction[2] = -1.0/sqrt(5.0);

    // New distance
    distance = 2.0;

    tracker->tally_crossings_implicit_abs(
	position, direction, distance, ew, sigma, tally);

    if (!soft_equiv(tally.get_outward_tally(2)[1], exp(-2.0) ) ) ITFAILS;
    if (!soft_equiv(tally.get_outward_tally(3)[1], exp(-1.0-sqrt(5.0)))) ITFAILS;
    if (!soft_equiv(tally.get_inward_tally(3) [1], exp(3.0-2*sqrt(5.0)))) ITFAILS;
    
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
