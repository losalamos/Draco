//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstExtrinsic_Tracker.cc
 * \author Mike Buksas
 * \date   Mon Jul 21 10:10:06 2003
 * \brief  unit test for Extrinsic_Surface_Tracker and
 *         Extrinsic_Tracker_Builder 
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "imc_test.hh"
#include "mc/RZWedge_Mesh.hh"
#include "mc/RZWedge_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Rep_Topology.hh"
#include "mc/Global_Mesh_Data.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "IMC_Test.hh"
#include "../Extrinsic_Surface_Tracker.hh"
#include "../Extrinsic_Tracker_Builder.hh"
#include "../Azimuthal_Mesh.hh"
#include "../Surface_Sub_Tally.hh"
#include "../Surface_Tracking_Interface.hh"

using namespace std;
using namespace rtt_mc;
using namespace rtt_imc;
using namespace rtt_imc_test;

using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef rtt_mc::RZWedge_Mesh RZ;

//---------------------------------------------------------------------------//
// Test interface implemetation
//---------------------------------------------------------------------------//

struct Surface_Tracking_Tester : public Surface_Tracking_Interface
{
    double small_radius, large_radius;

    vector<double> bin_cosines;
    vector<Surface_Descriptor> descriptor;

    Surface_Tracking_Tester(double small, double large);

    int number_of_surfaces() const { return 3; }

    const vector<Surface_Descriptor>& get_surface_data() const 
    { 
	return descriptor;
    }

    const vector<double>& get_bin_cosines() const { return bin_cosines; }

    ~Surface_Tracking_Tester() { /* ... */ }
};

Surface_Tracking_Tester::Surface_Tracking_Tester(double small, double large)
    : small_radius(small), large_radius(large)
{
    descriptor.resize(3);

    descriptor[0].type = Surface_Descriptor::SPHERE;
    descriptor[0].data.resize(2);
    descriptor[0].data[0] = 0.0;
    descriptor[0].data[1] = small_radius;
    
    descriptor[1].type = Surface_Descriptor::SPHERE;
    descriptor[1].data.resize(2);
    descriptor[1].data[0] = 1.0;
    descriptor[1].data[1] = small_radius;
    
//     descriptor[2].type = Surface_Descriptor::SPHERE;
//     descriptor[2].data.resize(2);
//     descriptor[2].data[0] = 0.0;
//     descriptor[2].data[1] = 50.0;
    
    descriptor[2].type = Surface_Descriptor::SPHERE;
    descriptor[2].data.resize(2);
    descriptor[2].data[0] = 2.0;
    descriptor[2].data[1] = large_radius;

    double bin_data[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    bin_cosines.assign(bin_data, bin_data+5);
}

struct Surface_Tracking_Tester_Zero : public Surface_Tracking_Interface
{
    vector<double> bin_cosines;
    vector<Surface_Descriptor> descriptor;

    int number_of_surfaces() const { return 0; }

    const vector<Surface_Descriptor>& get_surface_data() const 
    { 
	return descriptor;
    }

    const vector<double>& get_bin_cosines() const { return bin_cosines; }
};

//---------------------------------------------------------------------------//
// BUILDERS
//---------------------------------------------------------------------------//

Surface_Sub_Tally make_surface_tally(const Surface_Tracking_Interface& interface)
{
    SP<Azimuthal_Mesh> mesh ( new Azimuthal_Mesh( interface ) );
    return Surface_Sub_Tally(mesh, interface);
}

//---------------------------------------------------------------------------//

SP<RZWedge_Mesh> build_mesh()
{
    // make a builder from the RZWedge input
    SP<Parser> parser(new Parser("RZWedge_Input"));
    RZWedge_Builder builder(parser);

    SP<RZWedge_Mesh> mesh = builder.build_Mesh();

    Ensure(mesh);

    return mesh;
}

//---------------------------------------------------------------------------//

template<class MT>
SP<Global_Mesh_Data<MT> > build_Mesh_Data(const MT &mesh)
{
    SP<Topology> topology(new Rep_Topology(mesh.num_cells()));
    
    SP<Global_Mesh_Data<MT> > mesh_data(
	new Global_Mesh_Data<MT>(topology, mesh));

    Ensure (mesh_data);
    return mesh_data;
}

//---------------------------------------------------------------------------//

SP<Extrinsic_Tracker_Builder<RZ> > build_tracker_builder(
    const RZWedge_Mesh&                   mesh,
    const Global_Mesh_Data<RZWedge_Mesh>& mesh_data,
    const Surface_Tracking_Interface&     interface)
{
    SP<Extrinsic_Tracker_Builder<RZ> > builder (
	new Extrinsic_Tracker_Builder<RZ>(mesh, mesh_data, interface) );

    Ensure (builder);

    return builder;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_RZWedge_Mesh_tracker()
{
    SP<RZWedge_Mesh>                    mesh      = build_mesh();
    SP<Global_Mesh_Data<RZWedge_Mesh> > mesh_data = build_Mesh_Data(*mesh);

    double x1 = mesh->get_high_x(1);
    double x3 = mesh->get_high_x(3);
    Surface_Tracking_Tester tester(x1,x3);
    
    SP<Extrinsic_Tracker_Builder<RZ> > builder = 
	build_tracker_builder(*mesh, *mesh_data, tester);

    // Test the builder:
    if (builder->get_global_surfaces() != 3) ITFAILS;
    if (builder->get_local_surfaces()  != 3) ITFAILS;
    if (builder->get_cell_status(1) != true) ITFAILS;
    if (builder->get_cell_status(2) != true) ITFAILS;
    if (builder->get_cell_status(3) != true) ITFAILS;
    if (builder->get_cell_status(4) != true) ITFAILS;
    if (builder->get_cell_status(5) != true) ITFAILS;
    if (builder->get_cell_status(6) != true) ITFAILS;

    SP<Extrinsic_Surface_Tracker> tracker = builder->build_tracker();

    // Test the tracker
    if (tracker->surface_in_cell(1) != true) ITFAILS;
    if (tracker->surface_in_cell(2) != true) ITFAILS;
    if (tracker->surface_in_cell(3) != true) ITFAILS;
    if (tracker->surface_in_cell(4) != true) ITFAILS;
    if (tracker->surface_in_cell(5) != true) ITFAILS;
    if (tracker->surface_in_cell(6) != true) ITFAILS;

    vector<double> position(3);
    position[0]  = x1*0.5; position[1]  = 0.0; position[2]  = 0.0;

    vector<double> direction(3);
    direction[0] = 0.0;    direction[1] = 0.0; direction[2] = 1.0;

    double distance = 3.0;
    double sigma = 1.0;
    double ew = 1.0;

    int cell = mesh->get_cell(position);  
    if (cell != 1) ITFAILS;

    Surface_Sub_Tally tally ( make_surface_tally(tester) );

    tracker->initialize_status(position, direction);

    if (tracker->get_inside(1) != true)  ITFAILS;
    if (tracker->get_inside(2) != false) ITFAILS;
    if (tracker->get_inside(3) != true)  ITFAILS;

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

void test_OS_Mesh_tracker()
{
    // build a parser, mesh, and interface
    SP<Parser>      parser(new Parser("OS_Input_B"));
    SP<OS_Builder>  mb(new OS_Builder(parser));
    SP<OS_Mesh>     mesh = mb->build_Mesh();

    // build the mesh data
    SP<Global_Mesh_Data<OS_Mesh> > mesh_data = build_Mesh_Data(*mesh);

    double x1 = 1.0;
    double x3 = 2.0;
    Surface_Tracking_Tester tester(x1,x3);

    bool caught = false;
    try
    {
	Extrinsic_Tracker_Builder<OS_Mesh> builder(*mesh, *mesh_data, tester);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	ostringstream m;
	m << "Good, caught assertion, " << ass.what();
	PASSMSG(m.str());
	caught = true;
    }

    if (!caught)
    {
	FAILMSG("Failed to catch OS_Mesh assertion in builder.");
    }
}

//---------------------------------------------------------------------------//

void test_null_tracker()
{
    SP<RZWedge_Mesh>                    mesh      = build_mesh();
    SP<Global_Mesh_Data<RZWedge_Mesh> > mesh_data = build_Mesh_Data(*mesh);

    Surface_Tracking_Tester_Zero tester;
    
    SP<Extrinsic_Tracker_Builder<RZ> > builder = 
	build_tracker_builder(*mesh, *mesh_data, tester);

    SP<Extrinsic_Surface_Tracker> tracker = builder->build_tracker();

    // we shouldn't have a tracker because there are no global surfaces
    if (tracker) ITFAILS;
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
	test_RZWedge_Mesh_tracker();
	test_OS_Mesh_tracker();
	test_null_tracker();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstExtrinsic_Tracker, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_imc_test::passed) 
    {
        std::cout << "**** tstExtrinsic_Tracker Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstExtrinsic_Tracker." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstExtrinsic_Tracker.cc
//---------------------------------------------------------------------------//
