//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/MC_Test.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 11 15:22:12 2000
 * \brief  Components used for MC package testing
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../AMR_Layout.hh"
#include "../XYZCoord_sys.hh"
#include "../Constants.hh"
#include "c4/global.hh"
#include <iostream>
#include <string>
#include <cmath>

namespace rtt_mc_test
{

using rtt_mc::global::pi;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::Coord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Topology;
using rtt_mc::OS_Mesh;
using rtt_mc::AMR_Layout;
using rtt_dsxx::SP;
using std::vector;

using std::sqrt;
using std::cos;
using std::sin;

//===========================================================================//
// Topology class tests
//===========================================================================//
// Full replication topology test --> uses base class reference to derived
// class type

bool topology_replication_test(SP<OS_Mesh> mesh, const Topology &top)
{
    // passing condition
    bool p = true;

    // test num_cells
    if (top.num_cells() != mesh->num_cells()) p = FAILURE;
    if (top.num_cells(C4::node()) != mesh->num_cells()) p = FAILURE;
    
    // test num procs
    for (int cell = 1; cell <= top.num_cells(); cell++)
	if (top.num_procs(cell) != C4::nodes()) p = FAILURE;

    // verify the parallel scheme
    if (top.get_parallel_scheme() != "replication") p = FAILURE;

    // test global_cell functions
    for (int lc = 1; lc <= top.num_cells(C4::node()); lc++)
	if (top.global_cell(lc) != lc) p = FAILURE;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int lc = 1; lc <= top.num_cells(proc); lc++)
	    if (top.global_cell(lc, proc) != lc) p = FAILURE;

    // test local cell functions
    for (int gc = 1; gc <= top.num_cells(); gc++)
	if (top.local_cell(gc) != gc) p = FAILURE;

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int gc = 1; gc <= top.num_cells(); gc++)
	    if (top.local_cell(gc, proc) != gc) p = FAILURE;

    // test boundary cell functions
    for (int proc = 0; proc < C4::nodes(); proc++)
    {
	if (top.get_boundary_cells(proc) != 0) p = FAILURE; 
	for (int gc = 1; gc <= top.num_cells(); gc++)
	    if (top.global_to_boundary(gc, proc) != 0) p = FAILURE;
    }

    for (int proc = 0; proc < C4::nodes(); proc++)
	for (int bc = 1; bc <= top.get_boundary_cells(proc); bc++)
	{ 
	    // shouldn't get here because there are no boundary
	    // cells
	    p = FAILURE; 
	}

    // test get_cells function
    {
	vector<int> cell_list(mesh->num_cells());
	for (int i = 1; i <= cell_list.size(); i++)
	    cell_list[i-1] = i;

	// compare
	for (int np = 0; np < C4::nodes(); np++)
	{
	    vector<int> got_cells = top.get_cells(np);
	    if (got_cells != cell_list) p = FAILURE;
	}
    }

    // test get_procs function
    {
	vector<int> proc_list(C4::nodes());
	for (int i = 0; i < C4::nodes(); i++)
	    proc_list[i] = i;

	// compare
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    vector<int> got_procs = top.get_procs(cell);
	    if (got_procs != proc_list) p = FAILURE;
	}
    }

    // return passing condition
    return p;
}

//---------------------------------------------------------------------------//
// Full DD topology test --> uses base class reference to derived
// class type

bool topology_DD_test(SP<OS_Mesh> mesh, const Topology &top)
{
    // passing condition
    bool p = true;

    // verify the parallel scheme
    if (top.get_parallel_scheme() != "DD") p = FAILURE;

    // check the number of cells, this is a 2 processor/3 cell mesh
    // decomposition
    if (top.num_cells() != mesh->num_cells()) p = FAILURE;
    if (top.num_cells(C4::node()) != 3) p = FAILURE;

    // check the number of procs--each cell should be on one processor
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	if (top.num_procs(cell) != 1) p = FAILURE;

    // check global cell indices
    if (C4::node() == 0)
    {
	if (top.global_cell(1) != 1) p = FAILURE;
	if (top.global_cell(2) != 2) p = FAILURE;
	if (top.global_cell(3) != 3) p = FAILURE;
    }

    if (C4::node() == 1)
    {
	if (top.global_cell(1) != 4) p = FAILURE;
	if (top.global_cell(2) != 5) p = FAILURE;
	if (top.global_cell(3) != 6) p = FAILURE;
    }

    for (int np = 0; np < C4::nodes(); np++)
    {
	int offset = np * 3;
	for (int lc = 1; lc <= 3; lc++)
	    if (top.global_cell(lc, np) != lc + offset) p = FAILURE;
    }
	

    // check local cell indices
    if (C4::node() == 0)
    {
	if (top.local_cell(1) != 1) p = FAILURE;
	if (top.local_cell(2) != 2) p = FAILURE;
	if (top.local_cell(3) != 3) p = FAILURE;
	if (top.local_cell(4) != 0) p = FAILURE;
	if (top.local_cell(5) != 0) p = FAILURE;
	if (top.local_cell(6) != 0) p = FAILURE;
    }

    if (C4::node() == 1)
    {
	if (top.local_cell(1) != 0) p = FAILURE;
	if (top.local_cell(2) != 0) p = FAILURE;
	if (top.local_cell(3) != 0) p = FAILURE;
	if (top.local_cell(4) != 1) p = FAILURE;
	if (top.local_cell(5) != 2) p = FAILURE;
	if (top.local_cell(6) != 3) p = FAILURE;
    }

    {
	if (top.local_cell(1, 0) != 1) p = FAILURE;
	if (top.local_cell(2, 0) != 2) p = FAILURE;
	if (top.local_cell(3, 0) != 3) p = FAILURE;
	if (top.local_cell(4, 0) != 0) p = FAILURE;
	if (top.local_cell(5, 0) != 0) p = FAILURE;
	if (top.local_cell(6, 0) != 0) p = FAILURE;
	if (top.local_cell(1, 1) != 0) p = FAILURE;
	if (top.local_cell(2, 1) != 0) p = FAILURE;
	if (top.local_cell(3, 1) != 0) p = FAILURE;
	if (top.local_cell(4, 1) != 1) p = FAILURE;
	if (top.local_cell(5, 1) != 2) p = FAILURE;
	if (top.local_cell(6, 1) != 3) p = FAILURE;
    }

    // check boundary cells
    if (top.get_boundary_cells(0) != 3) p = FAILURE;
    if (top.get_boundary_cells(1) != 3) p = FAILURE;
    {
	if (top.global_to_boundary(1, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(2, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(3, 0) != 0) p = FAILURE;
	if (top.global_to_boundary(4, 0) != 1) p = FAILURE;
	if (top.global_to_boundary(5, 0) != 2) p = FAILURE;
	if (top.global_to_boundary(6, 0) != 3) p = FAILURE;
	if (top.global_to_boundary(1, 1) != 1) p = FAILURE;
	if (top.global_to_boundary(2, 1) != 2) p = FAILURE;
	if (top.global_to_boundary(3, 1) != 3) p = FAILURE;
	if (top.global_to_boundary(4, 1) != 0) p = FAILURE;
	if (top.global_to_boundary(5, 1) != 0) p = FAILURE;
	if (top.global_to_boundary(6, 1) != 0) p = FAILURE;
    }
    
    for (int np = 0; np < C4::nodes(); np++)
    {
	int offset;
	if (np == 0) offset = 3;
	if (np == 1) offset = 0;
	for (int bc = 1; bc <= top.get_boundary_cells(np); bc++)
	    if (top.boundary_to_global(bc, np) != bc + offset) p = FAILURE;
    }

    if (top.boundary_to_global(1, 0) != 4) p = FAILURE;
    if (top.boundary_to_global(2, 0) != 5) p = FAILURE;
    if (top.boundary_to_global(3, 0) != 6) p = FAILURE;
    if (top.boundary_to_global(1, 1) != 1) p = FAILURE;
    if (top.boundary_to_global(2, 1) != 2) p = FAILURE;
    if (top.boundary_to_global(3, 1) != 3) p = FAILURE;

    // check cell lists on processor
    vector<int> cell_list(3);
    if (C4::node() == 0) 
    {
	cell_list[0] = 1;
	cell_list[1] = 2;
	cell_list[2] = 3;
	if (top.get_cells(0) != cell_list) p = FAILURE;
    }
    
    if (C4::node() == 1) 
    {
	cell_list[0] = 4;
	cell_list[1] = 5;
	cell_list[2] = 6;
	if (top.get_cells(1) != cell_list) p = FAILURE;
    }

    // check proc lists for global cells
    vector<int> proc_list(1);
    proc_list[0] = 0;
    if (top.get_procs(1) != proc_list) p = FAILURE;
    if (top.get_procs(2) != proc_list) p = FAILURE;
    if (top.get_procs(3) != proc_list) p = FAILURE;
    proc_list[0] = 1;
    if (top.get_procs(4) != proc_list) p = FAILURE;
    if (top.get_procs(5) != proc_list) p = FAILURE;
    if (top.get_procs(6) != proc_list) p = FAILURE;

    // return passing condition
    return p;
}

//===========================================================================//
// BUILD AMR RZWEDGE_MESH
//===========================================================================//

SP<RZWedge_Mesh> make_RZWedge_Mesh_AMR(double phi)
{
    double rphi = phi * pi / 180.0;
    
    // make coord
    SP<Coord_sys> coord(new XYZCoord_sys());

    // make layout (single level mesh)
    AMR_Layout layout(12, 6);
    {
	layout.set_size(1, 6, 2);
	layout(1, 1, 1) = 1;
	layout(1, 2, 1) = 7;
	layout(1, 3, 1) = 1;
	layout(1, 4, 1) = 1;
	layout(1, 5, 1) = 1;
	layout(1, 6, 1) = 2;
	layout(1, 6, 2) = 4;
    
	layout(2, 1, 1) = 2;
	layout(2, 2, 1) = 4;
	layout(2, 3, 1) = 2;
	layout(2, 4, 1) = 2;
	layout(2, 5, 1) = 1;
	layout(2, 6, 1) = 3;
    
	layout(3, 1, 1) = 3;
	layout(3, 2, 1) = 5;
	layout(3, 3, 1) = 3;
	layout(3, 4, 1) = 3;
	layout(3, 5, 1) = 2;
	layout(3, 6, 1) = 6;
    
	layout(4, 1, 1) = 2;
	layout(4, 2, 1) = 8;
	layout(4, 3, 1) = 4;
	layout(4, 4, 1) = 4;
	layout(4, 5, 1) = 1;
	layout(4, 6, 1) = 5;
    
	layout(5, 1, 1) = 3;
	layout(5, 2, 1) = 8;
	layout(5, 3, 1) = 5;
	layout(5, 4, 1) = 5;
	layout(5, 5, 1) = 4;
	layout(5, 6, 1) = 6;
    
	layout.set_size(6, 2, 2);
	layout.set_size(6, 5, 2);
	layout(6, 1, 1) = 6;
	layout(6, 2, 1) = 9;
	layout(6, 2, 2) = 10;
	layout(6, 3, 1) = 6;
	layout(6, 4, 1) = 6;
	layout(6, 5, 1) = 3;
	layout(6, 5, 2) = 5;
	layout(6, 6, 1) = 0;
    
	layout(7, 1, 1) = 1;
	layout(7, 2, 1) = 0;
	layout(7, 3, 1) = 7;
	layout(7, 4, 1) = 7;
	layout(7, 5, 1) = 7;
	layout(7, 6, 1) = 8;
    
	layout.set_size(8, 1, 2);
	layout.set_size(8, 6, 2);
	layout(8, 1, 1) = 4;
	layout(8, 1, 2) = 5;
	layout(8, 2, 1) = 0;
	layout(8, 3, 1) = 8;
	layout(8, 4, 1) = 8;
	layout(8, 5, 1) = 7;
	layout(8, 6, 1) = 9;
	layout(8, 6, 2) = 11;
    
	layout(9, 1, 1) = 6;
	layout(9, 2, 1) = 11;
	layout(9, 3, 1) = 9;
	layout(9, 4, 1) = 9;
	layout(9, 5, 1) = 8;
	layout(9, 6, 1) = 10;
    
	layout(10, 1, 1) = 6;
	layout(10, 2, 1) = 12;
	layout(10, 3, 1) = 10;
	layout(10, 4, 1) = 10;
	layout(10, 5, 1) = 9;
	layout(10, 6, 1) = 0;
    
	layout(11, 1, 1) = 9;
	layout(11, 2, 1) = 0;
	layout(11, 3, 1) = 11;
	layout(11, 4, 1) = 11;
	layout(11, 5, 1) = 8;
	layout(11, 6, 1) = 12;
    
	layout(12, 1, 1) = 10;
	layout(12, 2, 1) = 0;
	layout(12, 3, 1) = 12;
	layout(12, 4, 1) = 12;
	layout(12, 5, 1) = 11;
	layout(12, 6, 1) = 0;
    }

    // make x-z extents for each cell
    vf_double xz(12, sf_double(4));
    vf_double rz(12, sf_double(4));
    {
	double rconv = sqrt(rphi/sin(rphi))*cos(rphi/2);

	rz[0][0] = 0/rconv;
	rz[0][1] = 1/rconv;
	rz[0][2] = 0;
	rz[0][3] = 1;
	
	rz[1][0] = 0/rconv;
	rz[1][1] = .5/rconv;
	rz[1][2] = 1;
	rz[1][3] = 1.5;
	
	rz[2][0] = 0/rconv;
	rz[2][1] = .5/rconv;
	rz[2][2] = 1.5;
	rz[2][3] = 2;
	
	rz[3][0] = .5/rconv;
	rz[3][1] = 1/rconv;
	rz[3][2] = 1;
	rz[3][3] = 1.5;
	
	rz[4][0] = .5/rconv;
	rz[4][1] = 1/rconv;
	rz[4][2] = 1.5;
	rz[4][3] = 2;
	
	rz[5][0] = 0/rconv;
	rz[5][1] = 1/rconv;
	rz[5][2] = 2;
	rz[5][3] = 3;
	
	rz[6][0] = 1/rconv;
	rz[6][1] = 2/rconv;
	rz[6][2] = 0;
	rz[6][3] = 1;
	
	rz[7][0] = 1/rconv;
	rz[7][1] = 2/rconv;
	rz[7][2] = 1;
	rz[7][3] = 2;
	
	rz[8][0] = 1/rconv;
	rz[8][1] = 1.5/rconv;
	rz[8][2] = 2;
	rz[8][3] = 2.5;
	
	rz[9][0] = 1/rconv;
	rz[9][1] = 1.5/rconv;
	rz[9][2] = 2.5;
	rz[9][3] = 3;
	
	rz[10][0] = 1.5/rconv;
	rz[10][1] = 2/rconv;
	rz[10][2] = 2;
	rz[10][3] = 2.5;
	
	rz[11][0] = 1.5/rconv;
	rz[11][1] = 2/rconv;
	rz[11][2] = 2.5;
	rz[11][3] = 3;

	for (int i = 0; i < layout.num_cells(); i++)
	{
	    // convert rz faces
	    for (int j = 0; j < 2; j++)
	    {
		xz[i][j]   = sqrt(rphi/sin(rphi)) * rz[i][j] * cos(rphi/2.0); 
		xz[i][j+2] = rz[i][j+2];
	    }
	}
    }

    // make the mesh
    SP<RZWedge_Mesh> mesh(new RZWedge_Mesh(coord, layout, xz, phi));
    Check (mesh->full_Mesh());

    return mesh;
}

} // end namespace rtt_mc_test

//---------------------------------------------------------------------------//
//                              end of MC_Test.cc
//---------------------------------------------------------------------------//
