//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstExtrinsic_Tracker.cc
 * \author Mike Buksas
 * \date   Mon Jul 21 10:10:06 2003
 * \brief  unit test for Extrinsic_Surface_Tracker and Extrinsic_Tracker_Builder
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
#include "mc/RZWedge_Mesh.hh"
#include "mc/RZWedge_Builder.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "IMC_Test.hh"
#include "../Extrinsic_Surface_Tracker.hh"
#include "../Extrinsic_Tracker_Builder.hh"
#include "../Azimuthal_Mesh.hh"
#include "../Surface_Sub_Tally.hh"

using namespace std;
using namespace rtt_mc;
using namespace rtt_imc;
using namespace rtt_imc_test;

using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// BUILDERS
//---------------------------------------------------------------------------//

Surface_Sub_Tally make_surface_tally()
{

    double c[3] = {-0.5, 0.0, 0.5};
    vector<double> cosines(c, c+3);

    SP<Azimuthal_Mesh> mesh ( new Azimuthal_Mesh( cosines ) );

    // The highest surface index number is 4:
    return Surface_Sub_Tally(mesh, 4);

}


SP<RZWedge_Mesh> build_mesh()
{

    // make a builder from the RZWedge input
    SP<Parser> parser(new Parser("RZWedge_Input"));
    RZWedge_Builder builder(parser);

    SP<RZWedge_Mesh> mesh = builder.build_Mesh();

    Ensure(mesh);

    return mesh;

}


SP<Extrinsic_Tracker_Builder> build_tracker_builder(const RZWedge_Mesh& mesh)
{

    double x1 = mesh.get_high_x(1);
    double x3 = mesh.get_high_x(3);

    SP<Extrinsic_Tracker_Builder> builder (
	new Extrinsic_Tracker_Builder(mesh) );

    builder->add_sphere(0.0, x1);    // Surface 1
    builder->add_sphere(1.0, x1);    // Surface 2
    builder->add_sphere(0.0, 50.0);  // Surface 3: Off the mesh
    builder->add_sphere(2.0, x3);    // Surface 4

    Ensure (builder);

    return builder;

}

SP<Extrinsic_Surface_Tracker> build_tracker(Extrinsic_Tracker_Builder& builder)
{

    SP<Extrinsic_Surface_Tracker> tracker = builder.build_tracker();

    Ensure(tracker);

    return tracker;

}





//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_tracker_builder()
{

    SP<RZWedge_Mesh> mesh = build_mesh();

    SP<Extrinsic_Tracker_Builder> builder = build_tracker_builder(*mesh);

    if (builder->get_global_surfaces() != 4) ITFAILS;

    if (builder->get_local_surfaces()  != 3) ITFAILS;

    if (builder->get_cell_status(1) != true)  ITFAILS;
    if (builder->get_cell_status(2) != true)  ITFAILS;
    if (builder->get_cell_status(3) != true)  ITFAILS;
    if (builder->get_cell_status(4) != true)  ITFAILS;
    if (builder->get_cell_status(5) != false) ITFAILS;
    if (builder->get_cell_status(6) != true)  ITFAILS;

}


void test_tracker()
{

    SP<RZWedge_Mesh> mesh = build_mesh();

    SP<Extrinsic_Tracker_Builder> builder = build_tracker_builder(*mesh);

    SP<Extrinsic_Surface_Tracker> tracker = builder->build_tracker();

    if (tracker->surface_in_cell(1) != true)  ITFAILS;
    if (tracker->surface_in_cell(2) != true)  ITFAILS;
    if (tracker->surface_in_cell(3) != true)  ITFAILS;
    if (tracker->surface_in_cell(4) != true)  ITFAILS;
    if (tracker->surface_in_cell(5) != false) ITFAILS;
    if (tracker->surface_in_cell(6) != true)  ITFAILS;

    double x1 = mesh->get_high_x(1);

    vector<double> position(3);
    position[0]  = x1*0.5; position[1]  = 0.0; position[2]  = 0.0;

    vector<double> direction(3);
    direction[0] = 0.0;    direction[1] = 0.0; direction[2] = 1.0;

    double distance = 3.0;
    double sigma = 1.0;
    double ew = 1.0;

    int cell = mesh->get_cell(position);  
    if (cell != 1) ITFAILS;

    Surface_Sub_Tally tally (make_surface_tally() );

    tracker->initialize_status(position, direction);

    if (tracker->get_inside(1) != true)  ITFAILS;
    if (tracker->get_inside(2) != false) ITFAILS;
    if (tracker->get_inside(4) != true)  ITFAILS;

    tracker->tally_crossings_implicit_abs(
	position, direction, cell, distance, ew, sigma, tally);

    if (!soft_equiv(tally.weight(1,true,4) , 
		    ew * exp(-sigma*x1*sqrt(3.0)/2.0) ) ) ITFAILS;

    if (!soft_equiv(tally.weight(2,false,4) ,
		    ew * exp(-sigma*(1.0-x1*sqrt(3.0)/2.0) ) ) ) ITFAILS;

    if (!soft_equiv(tally.weight(2,true,4) ,
		    ew * exp(-sigma*(1.0+x1*sqrt(3.0)/2.0) ) ) ) ITFAILS;

    

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
	test_tracker_builder();

	test_tracker();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstExtrinsic_Tracker, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" << std::endl;
    if (rtt_imc_test::passed) 
    {
        std::cout << "**** tstExtrinsic_Tracker Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstExtrinsic_Tracker." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstExtrinsic_Tracker.cc
//---------------------------------------------------------------------------//
