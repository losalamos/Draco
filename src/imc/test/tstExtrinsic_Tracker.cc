//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstExtrinsic_Tracker.cc
 * \author Thomas M. Evans
 * \date   Tue Jan  6 12:11:25 2004
 * \brief unit test for Extrinsic_Surface_Tracker and
 *        Extrinsic_Tracker_Builder
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
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "mc/RZWedge_Mesh.hh"
#include "mc/RZWedge_Builder.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/Rep_Topology.hh"
#include "mc/General_Topology.hh"
#include "mc/Global_Mesh_Data.hh"
#include "../Release.hh"
#include "../Extrinsic_Surface_Tracker.hh"
#include "../Extrinsic_Tracker_Builder.hh"
#include "../Azimuthal_Mesh.hh"
#include "../Surface_Sub_Tally.hh"
#include "../Surface_Tracking_Interface.hh"
#include "imc_test.hh"
#include "IMC_Test.hh"

using namespace std;

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
    
    descriptor[2].type = Surface_Descriptor::SPHERE;
    descriptor[2].data.resize(2);
    descriptor[2].data[0] = 2.0;
    descriptor[2].data[1] = large_radius;

    double bin_data[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    bin_cosines.assign(bin_data, bin_data+5);
}

//---------------------------------------------------------------------------//

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

struct Surface_Tracking_Tester_Sub_Mesh : public Surface_Tracking_Interface
{
    vector<double> bin_cosines;
    vector<Surface_Descriptor> descriptor;

    Surface_Tracking_Tester_Sub_Mesh();

    int number_of_surfaces() const { return 1; }

    const vector<Surface_Descriptor>& get_surface_data() const 
    { 
	return descriptor;
    }

    const vector<double>& get_bin_cosines() const { return bin_cosines; }

    ~Surface_Tracking_Tester_Sub_Mesh() { /* ... */ }
};

Surface_Tracking_Tester_Sub_Mesh::Surface_Tracking_Tester_Sub_Mesh()
    : descriptor(1)
{
    // need the sphere to reside entirely in cell 2 
    descriptor[0].type = Surface_Descriptor::SPHERE;
    descriptor[0].data.resize(2);
    descriptor[0].data[0] = 3.0;
    descriptor[0].data[1] = 0.5;

    double bin_data[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    bin_cosines.assign(bin_data, bin_data+5);
}

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

SP<RZWedge_Mesh> build_mesh_for_DD()
{
    Require (rtt_c4::nodes() == 2);

    // make a builder from the RZWedge input
    SP<Parser> parser(new Parser("RZWedge_Input_2_Cell"));
    RZWedge_Builder builder(parser);

    SP<RZWedge_Mesh> mesh = builder.build_Mesh();

    if (mesh->num_cells() != 2) ITFAILS;

    Ensure(mesh);

    return mesh;
}

//---------------------------------------------------------------------------//

SP<Topology> build_DD_topology()
{
    Require (rtt_c4::nodes() == 2);

    vector<vector<int> > cpp(2);
    vector<vector<int> > ppc(2);
    vector<vector<int> > bc(2);

    cpp[0].resize(1);
    cpp[1].resize(1);

    ppc[0].resize(1);
    ppc[1].resize(1);

    bc[0].resize(1);
    bc[1].resize(1);

    cpp[0][0] = 1;
    cpp[1][0] = 2;

    ppc[0][0] = 0;
    ppc[1][0] = 1;
    
    bc[0][0] = 2;
    bc[1][0] = 1;

    SP<Topology> top(new General_Topology(cpp, ppc, bc, "DD"));
    
    if (top->num_cells() != 2)  ITFAILS;
    if (top->num_cells(0) != 1) ITFAILS;
    if (top->num_cells(1) != 1) ITFAILS;

    if (top->global_cell(1, 0) != 1) ITFAILS;
    if (top->global_cell(1, 1) != 2) ITFAILS;

    if (top->boundary_to_global(1, 0) != 2) ITFAILS;
    if (top->boundary_to_global(1, 1) != 1) ITFAILS;

    return top;
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

void test_DD_tracker()
{
    if (rtt_c4::nodes() != 2) 
	return;

    SP<RZWedge_Mesh> mesh = build_mesh_for_DD();

    // put cell 1 on processor 0 and cell 2 on processor 1
    vector<int> map(2);
    if (rtt_c4::node() == 0)
    {
	map[0] = 1;
	map[1] = 0;
	SP<RZWedge_Mesh::Pack> pack = mesh->pack(map);
	mesh = pack->unpack();

	if (mesh->get_low_z(1) != 0.0)  ITFAILS;
	if (mesh->get_high_z(1) != 2.0) ITFAILS;
    }
    else if (rtt_c4::node() == 1)
    {
	map[0] = 0;
	map[1] = 1;
	SP<RZWedge_Mesh::Pack> pack = mesh->pack(map);
	mesh = pack->unpack();

	if (mesh->get_low_z(1) != 2.0)  ITFAILS;
	if (mesh->get_high_z(1) != 4.0) ITFAILS;
    }

    if (mesh->num_cells() != 1) ITFAILS;

    // make DD topology
    SP<Topology> dd = build_DD_topology();

    // make global mesh data
    Global_Mesh_Data<RZWedge_Mesh> mesh_data(dd, *mesh);

    // make interface that has 1 surface in global cell 2
    Surface_Tracking_Tester_Sub_Mesh interface;

    // make extrinsic tracker and builder
    Extrinsic_Tracker_Builder<RZWedge_Mesh> builder(*mesh, mesh_data, 
						    interface);

    // build the tracker
    SP<Extrinsic_Surface_Tracker> tracker = builder.build_tracker();

    if (!tracker)                                                   ITFAILS;
    if (tracker->get_num_global_surfaces() != 1)                    ITFAILS;
    if (!soft_equiv(tracker->get_surface_area(1), .0872665, 1.e-6)) ITFAILS;

    // now some on-proc checks
    if (rtt_c4::node() == 0)
    {
	if (builder.get_local_surfaces() != 0)      ITFAILS;
	if (builder.get_global_surfaces() != 1)     ITFAILS;

	if (tracker->surface_in_cell(1))            ITFAILS;

	if (tracker->surface_in_tracker(1))         ITFAILS;

	if (tracker->get_num_local_surfaces() != 0) ITFAILS;
    }

    if (rtt_c4::node() == 1)
    {
	if (builder.get_local_surfaces() != 1)      ITFAILS;
	if (builder.get_global_surfaces() != 1)     ITFAILS;

	if (!tracker->surface_in_cell(1))           ITFAILS;

	if (!tracker->surface_in_tracker(1))        ITFAILS;

	if (tracker->get_num_local_surfaces() != 1) ITFAILS;
    }

    if (rtt_imc_test::passed)
    {
	ostringstream m;
	m << "DD tracker tests ok on processor " << rtt_c4::node();
	PASSMSG(m.str());
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
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " 
		     << rtt_imc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	if (rtt_c4::node() == 0)
	{
	    test_RZWedge_Mesh_tracker();
	    test_OS_Mesh_tracker();
	    test_null_tracker();
	}

	// DD test
	test_DD_tracker();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstExtrinsic_Tracker, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_imc_test::passed) 
	{
	    cout << "**** tstExtrinsic_Tracker Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstExtrinsic_Tracker on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstExtrinsic_Tracker.cc
//---------------------------------------------------------------------------//
