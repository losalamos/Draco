//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstRZWedge_Mesh.cc
 * \author Todd J. Urbatsch
 * \date   Wed Apr 12 12:49:14 2000
 * \brief  Test the RZWedge_Mesh class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_tstRZWedge_Mesh_cc__
#define __test_tstRZWedge_Mesh_cc__
 
#include "MC_Test.hh"
#include "../RZWedge_Mesh.hh"
#include "../XYZCoord_sys.hh"
#include "../Layout.hh"
#include "../Release.hh"
#include "../Math.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "viz/Ensight_Translator.hh"
#include <iostream>
#include <vector>
#include "rng/Rnd_Control.hh"

using namespace std;

using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_dsxx::SP;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::global::soft_equiv;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void simple_one_cell_RZWedge()
{
    // >>> one-cell problem <<<
    int ncells = 1;

    // >>> build an XYZ coordinate system class <<<
    rtt_dsxx::SP<XYZCoord_sys> coord(new XYZCoord_sys());

    // >>> build a layout of size ncells <<<
    Layout layout;
    layout.set_size(ncells);

    // set number of faces per cell (6)
    for (int cell = 1; cell <= ncells; cell++)
	layout.set_size(cell, 6);

    // set vacuum boundary condition on external boundaries; 
    // reflecting on inherently reflecting boundaries (lox (radial axis),
    // loy, hiy)
    for (int cell = 1; cell <= ncells; cell++)
    {
	layout(cell, 1) = cell;
	layout(cell, 2) = 0; 
	layout(cell, 3) = cell;
	layout(cell, 4) = cell;
	layout(cell, 5) = 0;
	layout(cell, 6) = 0;
    }

    // >>> set XZ cell extents <<<
    vector<vector<double> > cell_xz_extents;
    cell_xz_extents.resize(ncells);
    for (int cell = 1; cell <= ncells; cell++)
	cell_xz_extents[cell-1].resize(4);

    for (int cell = 1; cell <= ncells; cell++)
    {
	cell_xz_extents[cell-1][0] = 0.0; // low x
	cell_xz_extents[cell-1][1] = 1.0; // high x
	cell_xz_extents[cell-1][2] = 0.0; // low z
	cell_xz_extents[cell-1][3] = 1.0; // high z
    }

    // set unfolding angle in degrees
    double theta_degrees = 90.0;
    bool submesh = false;

    // build a mesh object
    SP<RZWedge_Mesh> mesh(new RZWedge_Mesh(coord, layout, cell_xz_extents, 
    					   theta_degrees, submesh));

    // check that the mesh returns proper coordinate system information
    if (mesh->get_Coord().get_dim() != 3)    ITFAILS;
    if (mesh->get_SPCoord()->get_dim() != 3) ITFAILS;
    if (mesh->get_Coord().get_Coord()    != std::string("xyz")) ITFAILS;
    if (mesh->get_SPCoord()->get_Coord() != std::string("xyz")) ITFAILS;

    // check that mesh returns the proper number of cells
    if (mesh->num_cells() != ncells) ITFAILS;

    // check that x- and z-extents are correct
    if (mesh->get_low_x(1)  != 0.0)  ITFAILS;
    if (mesh->get_high_x(1) != 1.0)  ITFAILS;
    if (mesh->get_low_z(1)  != 0.0)  ITFAILS;
    if (mesh->get_high_z(1) != 1.0)  ITFAILS;

    // check that cell midpoints are correct
    if (mesh->get_x_midpoint(1) != 0.5) ITFAILS;
    if (mesh->get_y_midpoint(1) != 0.0) ITFAILS;
    if (mesh->get_z_midpoint(1) != 0.5) ITFAILS;

    // check that graphics cell type is correct
    vector<int> cell_type = mesh->get_cell_types();
    if (cell_type.size() != 1)                           ITFAILS;
    if (cell_type[0] != rtt_viz::eight_node_hexahedron ) ITFAILS;

    // check that mesh returns the proper point coordinates
    vector<vector<double> > point_coords = mesh->get_point_coord();
    if (point_coords.size() != 8) ITFAILS;
    if (point_coords[0].size() != 3) ITFAILS;
    if (point_coords[1].size() != 3) ITFAILS;
    if (point_coords[2].size() != 3) ITFAILS;
    if (point_coords[3].size() != 3) ITFAILS;
    if (point_coords[4].size() != 3) ITFAILS;
    if (point_coords[5].size() != 3) ITFAILS;
    if (point_coords[6].size() != 3) ITFAILS;
    if (point_coords[7].size() != 3) ITFAILS;

    if (!soft_equiv(point_coords[0][0], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[0][1], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[0][2], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[1][0], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[1][1],-1.0)) ITFAILS;
    if (!soft_equiv(point_coords[1][2], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[2][0], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[2][1], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[2][2], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[3][0], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[3][1], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[3][2], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[4][0], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[4][1], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[4][2], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[5][0], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[5][1],-1.0)) ITFAILS;
    if (!soft_equiv(point_coords[5][2], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[6][0], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[6][1], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[6][2], 1.0)) ITFAILS;
    if (!soft_equiv(point_coords[7][0], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[7][1], 0.0)) ITFAILS;
    if (!soft_equiv(point_coords[7][2], 1.0)) ITFAILS;

    // check that mesh returns the correct "next cell"
    if (mesh->next_cell(1, 1) != 1) ITFAILS;
    if (mesh->next_cell(1, 2) != 0) ITFAILS;
    if (mesh->next_cell(1, 3) != 1) ITFAILS;
    if (mesh->next_cell(1, 4) != 1) ITFAILS;
    if (mesh->next_cell(1, 5) != 0) ITFAILS;
    if (mesh->next_cell(1, 6) != 0) ITFAILS;

    // check that mesh can locate a cell given a position
    vector<double> position(3, 0.0);
    position[0] = 0.5;
    position[2] = 0.5;
    if (mesh->get_cell(position) != 1)  ITFAILS;

    position[0] = 100.0;
    if (mesh->get_cell(position) != -1) ITFAILS;

    // check that the face normals are correct
    vector<double> normal(3, 0.0);
    double sincos_45 = 1.0/std::sqrt(2.0);
    normal = mesh->get_normal(1,1);
    if (!soft_equiv(normal[0],-1.0)) ITFAILS;
    if (!soft_equiv(normal[1], 0.0)) ITFAILS;
    if (!soft_equiv(normal[2], 0.0)) ITFAILS;
    normal = mesh->get_normal(1,2);
    if (!soft_equiv(normal[0], 1.0)) ITFAILS;
    if (!soft_equiv(normal[1], 0.0)) ITFAILS;
    if (!soft_equiv(normal[2], 0.0)) ITFAILS;
    normal = mesh->get_normal(1,3);
    if (!soft_equiv(normal[0], -sincos_45)) ITFAILS;
    if (!soft_equiv(normal[1], -sincos_45)) ITFAILS;
    if (!soft_equiv(normal[2],        0.0)) ITFAILS;
    normal = mesh->get_normal(1,4);
    if (!soft_equiv(normal[0], -sincos_45)) ITFAILS;
    if (!soft_equiv(normal[1],  sincos_45)) ITFAILS;
    if (!soft_equiv(normal[2],        0.0)) ITFAILS;
    normal = mesh->get_normal(1,5);
    if (!soft_equiv(normal[0], 0.0)) ITFAILS;
    if (!soft_equiv(normal[1], 0.0)) ITFAILS;
    if (!soft_equiv(normal[2],-1.0)) ITFAILS;
    normal = mesh->get_normal(1,6);
    if (!soft_equiv(normal[0], 0.0)) ITFAILS;
    if (!soft_equiv(normal[1], 0.0)) ITFAILS;
    if (!soft_equiv(normal[2], 1.0)) ITFAILS;

    vector<double> normal_in(3, 0.0);
    for (int face = 1; face <= 6; face++)
    {
	normal = mesh->get_normal(1,face);
	normal_in = mesh->get_normal_in(1,face);
	for (int dim = 0; dim < 3; dim++)
	    if (!soft_equiv(normal[dim],-normal_in[dim])) ITFAILS;
    }

    // does the mesh return the correct cell volume?
    if (!soft_equiv(mesh->volume(1),1.0)) ITFAILS;

    // does the mesh return the correct face areas?
    if (!soft_equiv(mesh->face_area(1,1),0.0))            ITFAILS;
    if (!soft_equiv(mesh->face_area(1,2),2.0))            ITFAILS;
    if (!soft_equiv(mesh->face_area(1,3),std::sqrt(2.0))) ITFAILS;
    if (!soft_equiv(mesh->face_area(1,4),std::sqrt(2.0))) ITFAILS;
    if (!soft_equiv(mesh->face_area(1,5),1.0))            ITFAILS;
    if (!soft_equiv(mesh->face_area(1,6),1.0))            ITFAILS;


    // check that the mesh returns the correct vertices.
    // first, set reference vertices.
    vector<vector<double> > ref_vert;
    ref_vert.resize(3);
    for (int dim = 0; dim < 3; dim++)
    {
	ref_vert[dim].resize(8);
	for (int node = 0; node < 8; node++)
	    ref_vert[dim][node] = 0.0;
    }
    // node 1 (lox,loy,loz) - all zeros
    // node 2 (hix,loy,loz) 
    ref_vert[0][1] = 1.0;
    ref_vert[1][1] = -1.0;
    // node 3 (hix,hiy,loz) 
    ref_vert[0][2] = 1.0;
    ref_vert[1][2] = 1.0;
    // node 4 (lox,hiy,loz) - all zeros
    // node 5 (lox,loy,hiz)
    ref_vert[2][4] = 1.0;
    // node 6 (hix,loy,hiz) 
    ref_vert[0][5] = 1.0;
    ref_vert[1][5] = -1.0;
    ref_vert[2][5] = 1.0;
    // node 7 (hix,hiy,hiz) 
    ref_vert[0][6] = 1.0;
    ref_vert[1][6] = 1.0;
    ref_vert[2][6] = 1.0;
    // node 8 (lox,hiy,hiz) - all zeros
    ref_vert[2][7] = 1.0;

    vector<vector<int> > node_num(6,0);
    for (int face = 0; face < 6; face++)
	node_num[face].resize(4);
    // face 1 low x
    node_num[0][0] = 1;
    node_num[0][1] = 4;
    node_num[0][2] = 8;
    node_num[0][3] = 5;
    // face 2 high x
    node_num[1][0] = 2;
    node_num[1][1] = 3;
    node_num[1][2] = 7;
    node_num[1][3] = 6;
    // face 3 low y
    node_num[2][0] = 1;
    node_num[2][1] = 2;
    node_num[2][2] = 6;
    node_num[2][3] = 5;
    // face 4 high y
    node_num[3][0] = 4;
    node_num[3][1] = 3;
    node_num[3][2] = 7;
    node_num[3][3] = 8;
    // face 5 low z
    node_num[4][0] = 1;
    node_num[4][1] = 2;
    node_num[4][2] = 3;
    node_num[4][3] = 4;
    // face 6 high z
    node_num[5][0] = 5;
    node_num[5][1] = 6;
    node_num[5][2] = 7;
    node_num[5][3] = 8;


    vector<vector<double> > vertices;
    for (int face = 0; face < 6; face++)
    {
	vertices = mesh->get_vertices(1,face+1);

	if (vertices.size() != 3)    ITFAILS;  // 3 dimensions
	if (vertices[0].size() != 4) ITFAILS;  // 4 nodes per face
	if (vertices[1].size() != 4) ITFAILS;  // 4 nodes per face
	if (vertices[2].size() != 4) ITFAILS;  // 4 nodes per face

	for (int face_node = 0; face_node < 4; face_node++)
	    for (int dim = 0; dim < 3; dim++)
		if (!soft_equiv(vertices[dim][face_node],
				ref_vert[dim][node_num[face][face_node]-1]))
		    ITFAILS; 
    }

    vertices = mesh->get_vertices(1);

    if (vertices.size() != 3)    ITFAILS;  // 3 dimensions
    if (vertices[0].size() != 8) ITFAILS;  // 8 nodes per cell
    if (vertices[1].size() != 8) ITFAILS;  // 8 nodes per cell
    if (vertices[2].size() != 8) ITFAILS;  // 8 nodes per cell

    for (int cell_node = 0; cell_node < 4; cell_node++)
	for (int dim = 0; dim < 3; dim++)
	    if (!soft_equiv(vertices[dim][cell_node],
				ref_vert[dim][cell_node])) ITFAILS; 


    // check that the mesh returns the correct face number
    if (mesh->get_bndface("lox", 1) != 1) ITFAILS;
    if (mesh->get_bndface("hix", 1) != 2) ITFAILS;
    if (mesh->get_bndface("loy", 1) != 3) ITFAILS;
    if (mesh->get_bndface("hiy", 1) != 4) ITFAILS;
    if (mesh->get_bndface("loz", 1) != 5) ITFAILS;
    if (mesh->get_bndface("hiz", 1) != 6) ITFAILS;

    // check that, if this were a surface source cell, it would actually be
    // on a vacuume boundary
    vector<int> sur_cell(1,1);
    if ( mesh->check_defined_surcells("lox",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
    if ( mesh->check_defined_surcells("loy",sur_cell)) ITFAILS;
    if ( mesh->check_defined_surcells("hiy",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("loz",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("hiz",sur_cell)) ITFAILS;

    // check that the neighbors are correct
    vector<int> neighbors;
    neighbors = mesh->get_neighbors(1);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 1)     ITFAILS;
    if (neighbors[1] != 0)     ITFAILS;
    if (neighbors[2] != 1)     ITFAILS;
    if (neighbors[3] != 1)     ITFAILS;
    if (neighbors[4] != 0)     ITFAILS;
    if (neighbors[5] != 0)     ITFAILS;

    // make sure the mesh thinks it's not submesh
    if (!mesh->full_Mesh()) ITFAILS;

    // test the mesh's ability to sample positions.  This test has been hand
    // checked in a weak sense and, therefore, is more of a regression test.
    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    if (seed != 1234567) ITFAILS;

    vector<double> position_sampled;
    position_sampled = mesh->sample_pos(1,ran_object);
    if (position_sampled.size() != 3)                  ITFAILS;
    if (!soft_equiv(position_sampled[0],0.2515,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],0.1961,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],0.2216,0.001)) ITFAILS;

    position_sampled = mesh->sample_pos(1,ran_object);
    if (position_sampled.size() != 3)                   ITFAILS;
    if (!soft_equiv(position_sampled[0], 0.5232,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],-.02354,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2], 0.4629,0.001)) ITFAILS;

    vector<double> slope(3);
    slope[0] = -1.0;
    slope[1] =  0.0;
    slope[2] =  0.0;
    position_sampled = mesh->sample_pos(1, ran_object, slope, 1.0);
    if (position_sampled.size() != 3)                    ITFAILS;
    if (!soft_equiv(position_sampled[0],  0.2633,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],-.002651,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],  0.9373,0.001)) ITFAILS;

    
    position_sampled = mesh->sample_pos_on_face(1,2,ran_object);
    if (position_sampled.size() != 3)                   ITFAILS;
    if (!soft_equiv(position_sampled[0],   1.0,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],-.8054,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],0.3347,0.001)) ITFAILS;

    position_sampled = mesh->sample_pos_on_face(1,5,ran_object);
    if (position_sampled.size() != 3)                   ITFAILS;
    if (!soft_equiv(position_sampled[0], 0.4503,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],-.04137,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],    0.0,0.001)) ITFAILS;

    position_sampled = mesh->sample_pos_on_face(1,6,ran_object);
    if (position_sampled.size() != 3)                   ITFAILS;
    if (!soft_equiv(position_sampled[0], 0.2198,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1], 0.1769,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],    1.0,0.001)) ITFAILS;

    // check the == and != operations.
    // first, build another mesh object equivalent to old mesh object.
    // these two meshes, although identical, should not occupy the same memory
    SP<RZWedge_Mesh> other_mesh(new RZWedge_Mesh(coord, layout,
						 cell_xz_extents,  
						 theta_degrees, 
						 submesh));

    if (mesh == other_mesh) ITFAILS;

    RZWedge_Mesh &mesh_obj       = *mesh;
    RZWedge_Mesh &other_mesh_obj = *other_mesh;

    if (mesh_obj != other_mesh_obj)    ITFAILS;
    if (!(mesh_obj == other_mesh_obj)) ITFAILS;

    // check that the "other" mesh returns proper coordinate system information
    if (mesh->get_Coord().get_dim()    != 
	other_mesh->get_Coord().get_dim())	ITFAILS; 
    if (mesh->get_SPCoord()->get_dim() !=
	other_mesh->get_SPCoord()->get_dim())   ITFAILS;
    if (mesh->get_Coord().get_Coord()    != 
	other_mesh->get_Coord().get_Coord())    ITFAILS;
    if (mesh->get_SPCoord()->get_Coord() != 
	other_mesh->get_SPCoord()->get_Coord()) ITFAILS;


    // test the mesh's distance-to-boundary function
    vector<double> r(3, 0.0);
    vector<double> omega(3, 0.0);
    int intersecting_face;
    double db;
    
    r[2] = 0.5;
    omega[0] = 1.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db,1.0))    ITFAILS;
    if (intersecting_face != 2) ITFAILS;

    r[0] =  0.5;
    r[1] = -0.5;
    r[2] = 0.67823;
    omega[0] = 0.0;
    omega[1] = 1.0;
    omega[2] = 0.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db,1.0))    ITFAILS;
    if (intersecting_face != 4) ITFAILS;
    
    r[0] = 0.5;
    r[1] = 0.0;
    r[2] = 0.2398;
    omega[0] = 0.0;
    omega[1] = 0.0;
    omega[2] = 1.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db,0.7602)) ITFAILS;
    if (intersecting_face != 6) ITFAILS;

    r[0] = 0.5;
    r[1] = 0.0;
    r[2] = 0.7;
    omega[0] = -0.70710678119;
    omega[1] =  0.70710678119;
    omega[2] =  0.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db,0.35355339059, 1.0e-10)) ITFAILS;
    if (intersecting_face != 4)                 ITFAILS;
 
    r[0] = 0.5;
    r[1] = 0.0;
    r[2] = 0.49;
    omega[0] = -0.70710678119;
    omega[1] = -0.70710678119;
    omega[2] =  0.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db,0.35355339059, 1.0e-10)) ITFAILS;
    if (intersecting_face != 3)                 ITFAILS;

    r[0] = 0.5678;
    r[1] = 0.0;
    r[2] = 0.3;
    omega[0] = -1.0;
    omega[1] =  0.0;
    omega[2] =  0.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db, 0.5678, 1.0e-10))                    ITFAILS;
    if (!((intersecting_face == 1)||(intersecting_face == 3)
	 ||(intersecting_face == 4)))                        ITFAILS;

    r[0] =  1.0;
    r[1] = -0.5;
    r[2] =  0.1;
    omega[0] = -1.0;
    omega[1] =  0.0;
    omega[2] =  0.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db, 0.5, 1.0e-10)) ITFAILS;
    if (intersecting_face != 3)        ITFAILS;

    r[0] =  0.5;
    r[1] =  0.0;
    r[2] =  0.99;
    omega[0] =  0.0;
    omega[1] =  0.0;
    omega[2] = -1.0;
    db = mesh->get_db(r, omega, 1, intersecting_face);
    if (!soft_equiv(db, 0.99, 1.0e-10)) ITFAILS;
    if (intersecting_face != 5)         ITFAILS;

}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_mc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// simple one-celled test problem
	simple_one_cell_RZWedge();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing RZWedge_Mesh, " << ass.what() << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "****************************************" << endl;
    if (passed) 
    {
        cout << "**** RZWedge_Mesh Self Test: PASSED ****" << endl;
    }
    cout <<     "****************************************" << endl;
    cout << endl;

    cout << "Done testing RZWedge_Mesh." << endl;
    
    C4::Finalize();
}   

#endif                          // __test_tstRZWedge_Mesh_cc__

//---------------------------------------------------------------------------//
//                        end of test/tstRZWedge_Mesh.cc
//---------------------------------------------------------------------------//
