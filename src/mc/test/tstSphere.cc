//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSphere.cc
 * \author Mike Buksas
 * \date   Tue Jun 17 13:03:48 2003
 * \brief  Unit test for the Sphere implementation of the Surface interface
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include "../Release.hh"
#include "../Sphere.hh" 
#include "../Constants.hh"
#include "mc_test.hh"

using namespace std;
using namespace rtt_mc;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

SP<Sphere> make_sphere(double center, double raidus)
{
    return SP<Sphere> ( new Sphere(center, raidus) );
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_properties()
{

    const double center = 10.0;
    const double radius = 2.0;

    SP<Sphere> sphere ( make_sphere(center, radius) );

    if (!soft_equiv(sphere->surface_area(), 16.0*global::pi) ) ITFAILS;
    if (!soft_equiv(sphere->volume(), 32.0/3.0*global::pi) ) ITFAILS;

}

//---------------------------------------------------------------------------//

void test_distance_to()
{

    const double center = 0.0;
    const double raidus = 2.0;

    SP<Sphere> sphere( make_sphere(center, raidus) );

    vector<double> position(3);
    vector<double> direction(3);

    position[0]  = -3.0;    position[1]  = 0.0;    position[2]  = 0.0;
    direction[0] =  1.0;    direction[1] = 0.0;    direction[2] = 0.0;

    if (!soft_equiv(sphere->distance_to(position, direction), 1.0) ) ITFAILS;


    position[0]  =  0.0;    position[1]  = 0.0;    position[2]  = 0.0;
    direction[0] =  1.0;    direction[1] = 0.0;    direction[2] = 0.0;

    if (!soft_equiv(sphere->distance_to(position, direction), 2.0) ) ITFAILS;
    

    position[0]  =  3.0;    position[1]  = 0.0;    position[2]  = 0.0;
    direction[0] =  0.0;    direction[1] = 0.0;    direction[2] = 1.0;

    if ( !soft_equiv(sphere->distance_to(position, direction), 
		     global::huge) ) ITFAILS;


    position[0]  = -1.0;    position[1]  = 0.0;    position[2]  = 0.0;
    direction[0] = 1.0/std::sqrt(5.0);    
    direction[1] = 2.0/std::sqrt(5.0);    
    direction[2] = 0.0;

    if ( !soft_equiv(sphere->distance_to(position, direction), 
		     std::sqrt(5.0) ) ) ITFAILS;


}

//---------------------------------------------------------------------------//

void test_constrained_distance_to()
{

    const double center = 0.0;
    const double raidus = 2.0;

    SP<Sphere> sphere( make_sphere(center, raidus) );

    vector<double> position (3);
    vector<double> direction(3); 
    bool is_in;
    
    const double epsilon = 1.0e-7;
    double distance;

    position[0] = -2.0-epsilon;
    position[1] = 0.0;
    position[2] = 0.0;

    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;

    is_in = true;
    distance = sphere->distance_to(position, direction, is_in);
    if ( !soft_equiv(distance, 4.0, epsilon) ) ITFAILS;



    position[0] = -2.0+epsilon;
    distance = sphere->distance_to(position, direction, is_in);
    if ( !soft_equiv(distance, 4.0, epsilon) ) ITFAILS;



    is_in = false;
    try 
    {
	distance = sphere->distance_to(position, direction, is_in);
	ITFAILS;
    }
    catch (...) { }


    position[0] = -2.0-epsilon;
    distance = sphere->distance_to(position, direction, is_in);
    if ( !soft_equiv(distance, epsilon, epsilon) ) ITFAILS;



    position[0] = -2.0+epsilon;
    try 
    {
	distance = sphere->distance_to(position, direction, is_in);
	ITFAILS;
    }
    catch (...) { }

}

//---------------------------------------------------------------------------//

void test_is_inside()
{

    const double center = 0.0;
    const double radius = 2.0;

    SP<Sphere> sphere( make_sphere(center, radius) );

    vector<double> position (3);

    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;

    if (sphere->is_inside(position) != true) ITFAILS;

    position[1] = 3.0;

    if (sphere->is_inside(position) != false) ITFAILS;

    position[1] = 2.0;

    if (sphere->is_inside(position) != false) ITFAILS;

    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = center + radius;

    if (sphere->is_inside(position) != false) ITFAILS;

    vector<double> direction(3);
    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;

    if (sphere->is_inside(position, direction) != false) ITFAILS;

    direction[2] = -1.0;

    if (sphere->is_inside(position, direction) != true) ITFAILS;
    

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
	test_properties();
	
	test_distance_to();

	test_constrained_distance_to();

	test_is_inside();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstSphere, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_mc_test::passed) 
    {
        std::cout << "**** tstSphere Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstSphere." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSphere.cc
//---------------------------------------------------------------------------//
