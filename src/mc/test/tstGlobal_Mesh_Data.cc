//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstGlobal_Mesh_Data.cc
 * \author Thomas M. Evans
 * \date   Fri Dec  5 11:09:06 2003
 * \brief  test of Global_Mesh_Data class.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Soft_Equivalence.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../OS_Builder.hh"
#include "../OS_Mesh.hh"
#include "../RZWedge_Mesh.hh"
#include "../RZWedge_Builder.hh"
#include "../Sphyramid_Mesh.hh"
#include "../Sphyramid_Builder.hh"
#include "../Global_Mesh_Data.hh"
#include "../General_Topology.hh"
#include "../Rep_Topology.hh"
#include "../Release.hh"
#include "../Constants.hh"
#include "mc_test.hh"
#include "MC_Test.hh"

using namespace std;

using rtt_mc_test::Parser;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::RZWedge_Builder;
using rtt_mc::Sphyramid_Builder;
using rtt_mc::Sphyramid_Mesh;
using rtt_mc::General_Topology;
using rtt_mc::Rep_Topology;
using rtt_mc::Global_Mesh_Data;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

typedef vector<int>          sf_int;
typedef vector<vector<int> > vf_int;

//---------------------------------------------------------------------------//
// SERVICES
//---------------------------------------------------------------------------//

SP<Rep_Topology> build_rep_topology(int num_cells)
{
    SP<Rep_Topology> topology(new Rep_Topology(num_cells));
    return topology;
}

//---------------------------------------------------------------------------//
// for these tests we don't need a filled out topology

SP<General_Topology> build_general_topology()
{
    // for these
    // int nc  = 1;
    
    vf_int cpp(rtt_c4::nodes());
    vf_int ppc;
    vf_int bc(rtt_c4::nodes());
    
    SP<General_Topology> top(new General_Topology(cpp,ppc,bc,"DD"));
    return top;
}

//---------------------------------------------------------------------------//
// OS_Mesh builds

SP<OS_Mesh> build_OS_Mesh()
{
    SP<Parser>  parser(new Parser("OS_Input_3D"));
    OS_Builder  builder(parser);
    SP<OS_Mesh> mesh = builder.build_Mesh();
    if (mesh->num_cells() != 12) ITFAILS;
    return mesh;
}

//---------------------------------------------------------------------------//

SP<OS_Mesh> build_OS_Mesh_DD()
{
    Require (rtt_c4::nodes() == 4);
    
    SP<OS_Mesh> mesh;
    sf_int      c(12);
    if (rtt_c4::node() == 0)
    {
	mesh = build_OS_Mesh();

	// send cells 1-6 to processor 2
	{
	    for (int i = 0; i < 12; i++)
		c[i] = 0;
	    c[0] = 1;
	    c[1] = 2;
	    c[2] = 3;
	    c[3] = 4;
	    c[4] = 5;
	    c[5] = 6;
	    
	    SP<OS_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 1, 10);
	    rtt_c4::send(p->begin(), size, 1, 11);
	}

	// send cells 7 and 8 to processor 3
	{
	    for (int i = 0; i < 12; i++)
		c[i] = 0;
	    c[6] = 1;
	    c[7] = 2; 
	    
	    SP<OS_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 2, 10);
	    rtt_c4::send(p->begin(), size, 2, 11);
	}

	// send cells 9, 12 to processor 4
	{
	    for (int i = 0; i < 12; i++)
		c[i] = 0;
	    c[8]  = 1;
	    c[11] = 2; 
	    
	    SP<OS_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 3, 10);
	    rtt_c4::send(p->begin(), size, 3, 11);
	}

	// make processor 0 have cells 10, 11
	{
	    for (int i = 0; i < 12; i++)
		c[i] = 0;
	    c[9]  = 1; 
	    c[10] = 2;
	    
	    SP<OS_Mesh::Pack> p = mesh->pack(c);
	    mesh                = p->unpack();
	    if (mesh->num_cells() != 2) ITFAILS;

	    if (!soft_equiv(mesh->begin(1), -1.0)) ITFAILS;
	    if (!soft_equiv(mesh->end(1), 1.0))    ITFAILS;
	    if (!soft_equiv(mesh->begin(2), 1.0))  ITFAILS;
	    if (!soft_equiv(mesh->end(2), 3.0))    ITFAILS;
	    if (!soft_equiv(mesh->begin(3), 1.0))  ITFAILS;
	    if (!soft_equiv(mesh->end(3), 2.0))    ITFAILS;
	}
    }

    if (rtt_c4::node() == 1)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	OS_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 6) ITFAILS;
	
	if (!soft_equiv(mesh->begin(1), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->end(1), 2.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(2), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->end(2), 3.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(3), 0.0))  ITFAILS;
	if (!soft_equiv(mesh->end(3), 1.0))    ITFAILS;
    }

    if (rtt_c4::node() == 2)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	OS_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 2) ITFAILS;
	
	if (!soft_equiv(mesh->begin(1), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->end(1), 1.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(2), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->end(2), 1.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(3), 1.0))  ITFAILS;
	if (!soft_equiv(mesh->end(3), 2.0))    ITFAILS;
    }

    if (rtt_c4::node() == 3)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	OS_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 2) ITFAILS;
	
	if (!soft_equiv(mesh->begin(1), 1.0))  ITFAILS;
	if (!soft_equiv(mesh->end(1), 2.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(2), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->end(2), 3.0))    ITFAILS;
	if (!soft_equiv(mesh->begin(3), 1.0))  ITFAILS;
	if (!soft_equiv(mesh->end(3), 2.0))    ITFAILS;
    }

    return mesh;
}

//---------------------------------------------------------------------------//
// RZWedge_Mesh builds

SP<RZWedge_Mesh> build_RZWedge_Mesh()
{
    SP<Parser>       parser(new Parser("RZWedge_Input"));
    RZWedge_Builder  builder(parser);
    SP<RZWedge_Mesh> mesh = builder.build_Mesh();
    if (mesh->num_cells() != 6) ITFAILS;
    return mesh;   
}

//---------------------------------------------------------------------------//

SP<RZWedge_Mesh> build_RZWedge_Mesh_DD()
{
    Require (rtt_c4::nodes() == 4);

    SP<RZWedge_Mesh> mesh;
    sf_int           c(6);
    if (rtt_c4::node() == 0)
    {
	mesh = build_RZWedge_Mesh();

	// send cells 1-2 to processor 2
	{
	    for (int i = 0; i < 6; i++)
		c[i] = 0;
	    c[0] = 1;
	    c[1] = 2;
	    
	    SP<RZWedge_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 1, 10);
	    rtt_c4::send(p->begin(), size, 1, 11);
	}

	// send cells 3-4 to processor 3
	{
	    for (int i = 0; i < 6; i++)
		c[i] = 0;
	    c[2] = 1;
	    c[3] = 2;
	    
	    SP<RZWedge_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 2, 10);
	    rtt_c4::send(p->begin(), size, 2, 11);
	}

	// send cell 5 to processor 4
	{
	    for (int i = 0; i < 6; i++)
		c[i] = 0;
	    c[4] = 1;
	    
	    SP<RZWedge_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 3, 10);
	    rtt_c4::send(p->begin(), size, 3, 11);
	}

	// make processor 0 have cell 6
	{
	    for (int i = 0; i < 6; i++)
		c[i] = 0;
	    c[5]  = 1; 
	    
	    SP<RZWedge_Mesh::Pack> p = mesh->pack(c);
	    mesh                     = p->unpack();
	    if (mesh->num_cells() != 1) ITFAILS;

	    if (!soft_equiv(mesh->get_low_z(1), 1.0))  ITFAILS;
	    if (!soft_equiv(mesh->get_high_z(1), 3.0)) ITFAILS;
	}
    }

    if (rtt_c4::node() == 1)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	RZWedge_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 2) ITFAILS;
	
	if (!soft_equiv(mesh->get_low_x(1), 0.0))  ITFAILS;
	if (!soft_equiv(mesh->get_low_z(1), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->get_high_z(1), 1.0)) ITFAILS;
    }

    if (rtt_c4::node() == 2)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	RZWedge_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 2) ITFAILS;
	
	if (!soft_equiv(mesh->get_low_z(1), -1.0)) ITFAILS;
	if (!soft_equiv(mesh->get_high_z(1), 1.0)) ITFAILS;
    }

    if (rtt_c4::node() == 3)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	RZWedge_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 1) ITFAILS;
	
	if (!soft_equiv(mesh->get_low_z(1), 1.0))  ITFAILS;
	if (!soft_equiv(mesh->get_high_z(1), 3.0)) ITFAILS;
    }
    
    return mesh;
}

//---------------------------------------------------------------------------//
// Sphyramid_Mesh builds

SP<Sphyramid_Mesh> build_Sphyramid_Mesh()
{
    SP<Parser>         parser(new Parser("Sphyramid_Input"));
    Sphyramid_Builder  builder(parser);
    SP<Sphyramid_Mesh> mesh = builder.build_Mesh();
    if (mesh->num_cells() != 5) ITFAILS;
    return mesh;   
}

//---------------------------------------------------------------------------//

SP<Sphyramid_Mesh> build_Sphyramid_Mesh_DD()
{
    Require (rtt_c4::nodes() == 4);

    SP<Sphyramid_Mesh> mesh;
    sf_int             c(5);
    if (rtt_c4::node() == 0)
    {
	mesh = build_Sphyramid_Mesh();
	if (mesh->num_cells() != 5) ITFAILS;

	// send cells 1-2 to processor 2
	{
	    for (int i = 0; i < 5; i++)
		c[i] = 0;
	    c[0] = 1;
	    c[1] = 2;
	    
	    SP<Sphyramid_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 1, 10);
	    rtt_c4::send(p->begin(), size, 1, 11);
	}

	// send cell 3 to processor 3
	{
	    for (int i = 0; i < 5; i++)
		c[i] = 0;
	    c[2] = 1;
	    
	    SP<Sphyramid_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 2, 10);
	    rtt_c4::send(p->begin(), size, 2, 11);
	}

	// send cell 4 to processor 4
	{
	    for (int i = 0; i < 5; i++)
		c[i] = 0;
	    c[3] = 1;
	    
	    SP<Sphyramid_Mesh::Pack> p = mesh->pack(c);
	    int size = p->get_size();
	    rtt_c4::send(&size, 1, 3, 10);
	    rtt_c4::send(p->begin(), size, 3, 11);
	}

	// make processor 0 have cell 5
	{
	    for (int i = 0; i < 5; i++)
		c[i] = 0;
	    c[4]  = 1; 
	    
	    SP<Sphyramid_Mesh::Pack> p = mesh->pack(c);
	    mesh                       = p->unpack();
	    if (mesh->num_cells() != 1) ITFAILS;
	}
    }

    if (rtt_c4::node() == 1)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	Sphyramid_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 2) ITFAILS;
    }

    if (rtt_c4::node() == 2)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	Sphyramid_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 1) ITFAILS;
    }

    if (rtt_c4::node() == 3)
    {
	int   size;
	char *data;
	
	rtt_c4::receive(&size, 1, 0, 10);
	data = new char[size];
	rtt_c4::receive(data, size, 0, 11);

	Sphyramid_Mesh::Pack p(size, data);
	mesh = p.unpack();

	if (mesh->num_cells() != 1) ITFAILS;
    }
    
    return mesh;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_OS_Mesh_rep()
{
    SP<OS_Mesh>      mesh     = build_OS_Mesh();
    SP<Rep_Topology> topology = build_rep_topology(12);

    Global_Mesh_Data<OS_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 6) ITFAILS;

    // check x
    if (!soft_equiv(extents[0], -1.0)) ITFAILS;
    if (!soft_equiv(extents[1], 2.0))  ITFAILS;

    // check y
    if (!soft_equiv(extents[2], -1.0)) ITFAILS;
    if (!soft_equiv(extents[3], 3.0))  ITFAILS;

    // check z
    if (!soft_equiv(extents[4], 0.0))  ITFAILS;
    if (!soft_equiv(extents[5], 2.0))  ITFAILS;
    

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for OS_Mesh replication ok.");
}

//---------------------------------------------------------------------------//

void test_OS_Mesh_dd()
{
    if (rtt_c4::nodes() != 4)
	return;

    SP<OS_Mesh>          mesh     = build_OS_Mesh_DD();
    SP<General_Topology> topology = build_general_topology();

    Global_Mesh_Data<OS_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 6) ITFAILS;

    // check x
    if (!soft_equiv(extents[0], -1.0)) ITFAILS;
    if (!soft_equiv(extents[1], 2.0))  ITFAILS;

    // check y
    if (!soft_equiv(extents[2], -1.0)) ITFAILS;
    if (!soft_equiv(extents[3], 3.0))  ITFAILS;

    // check z
    if (!soft_equiv(extents[4], 0.0))  ITFAILS;
    if (!soft_equiv(extents[5], 2.0))  ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for OS_Mesh DD ok.");
}

//---------------------------------------------------------------------------//

void test_RZWedge_Mesh_rep()
{
    SP<RZWedge_Mesh> mesh     = build_RZWedge_Mesh();
    SP<Rep_Topology> topology = build_rep_topology(6);

    Global_Mesh_Data<RZWedge_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 4) ITFAILS;

    // high x
    double w  = sqrt(0.5 * rtt_mc::global::pi) * 3.0;
    double hx = w * sqrt(2.0) / 2.0;

    // check x
    if (!soft_equiv(extents[0], 0.0)) ITFAILS;
    if (!soft_equiv(extents[1], hx))  ITFAILS;

    // check z
    if (!soft_equiv(extents[2], -1.0)) ITFAILS;
    if (!soft_equiv(extents[3], 3.0))  ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for RZWedge_Mesh replication ok.");
}

//---------------------------------------------------------------------------//

void test_RZWedge_Mesh_dd()
{
    if (rtt_c4::nodes() != 4)
	return;

    SP<RZWedge_Mesh>     mesh     = build_RZWedge_Mesh_DD();
    SP<General_Topology> topology = build_general_topology();

    if (mesh->num_cells() == 6) ITFAILS;

    Global_Mesh_Data<RZWedge_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 4) ITFAILS;

    // high x
    double w  = sqrt(0.5 * rtt_mc::global::pi) * 3.0;
    double hx = w * sqrt(2.0) / 2.0;

    // check x
    if (!soft_equiv(extents[0], 0.0)) ITFAILS;
    if (!soft_equiv(extents[1], hx))  ITFAILS;

    // check z
    if (!soft_equiv(extents[2], -1.0)) ITFAILS;
    if (!soft_equiv(extents[3], 3.0))  ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for RZWedge_Mesh DD ok.");
}

//---------------------------------------------------------------------------//

void test_Sphyramid_Mesh_rep()
{
    SP<Sphyramid_Mesh> mesh     = build_Sphyramid_Mesh();
    SP<Rep_Topology>   topology = build_rep_topology(5);

    Global_Mesh_Data<Sphyramid_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 2) ITFAILS;

    // high x
    double rad = 45.0 * rtt_mc::global::pi / 180.0;
    double hx  = pow((2.0 * (1.0 - cos(rad)) / (tan(rad) * tan(rad))),
		     (1.0/3.0)) * 3.0;

    // check x
    if (!soft_equiv(extents[0], 0.0)) ITFAILS;
    if (!soft_equiv(extents[1], hx))  ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for Sphyramid_Mesh replication ok.");
}

//---------------------------------------------------------------------------//

void test_Sphyramid_Mesh_dd()
{
    if (rtt_c4::nodes() != 4)
	return;

    SP<Sphyramid_Mesh> mesh       = build_Sphyramid_Mesh_DD();
    SP<General_Topology> topology = build_general_topology();

    Global_Mesh_Data<Sphyramid_Mesh> gmd(topology, *mesh);

    // check spatial extents
    vector<double> extents = gmd.get_spatial_extents();
    if (extents.size() != 2) ITFAILS;

    // high x
    double rad = 45.0 * rtt_mc::global::pi / 180.0;
    double hx  = pow((2.0 * (1.0 - cos(rad)) / (tan(rad) * tan(rad))),
		     (1.0/3.0)) * 3.0;

    // check x
    if (!soft_equiv(extents[0], 0.0)) ITFAILS;
    if (!soft_equiv(extents[1], hx))  ITFAILS;

    if (rtt_mc_test::passed)
	PASSMSG("Global_Mesh_Data for Sphyramid_Mesh DD ok.");
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
		     << rtt_mc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_OS_Mesh_rep();
	test_OS_Mesh_dd();
	test_RZWedge_Mesh_rep();
	test_RZWedge_Mesh_dd();
	test_Sphyramid_Mesh_rep();
	test_Sphyramid_Mesh_dd();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstGlobal_Mesh_Data, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstGlobal_Mesh_Data Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstGlobal_Mesh_Data on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstGlobal_Mesh_Data.cc
//---------------------------------------------------------------------------//
