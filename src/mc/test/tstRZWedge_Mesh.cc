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
#include "../AMR_Layout.hh"
#include "../RZWedge_Builder.hh"
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

using rtt_mc_test::Parser_RZ;
using rtt_mc_test::make_RZWedge_Mesh_AMR;
using rtt_mc::XYZCoord_sys;
using rtt_mc::AMR_Layout;
using rtt_dsxx::SP;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::RZWedge_Builder;
using rtt_mc::global::soft_equiv;
using rtt_mc::global::pi;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;

using std::pow;

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
    AMR_Layout layout;
    layout.set_size(ncells);

    // set number of faces per cell (6)
    for (int cell = 1; cell <= ncells; cell++)
	layout.set_size(cell, 6);

    // set vacuum boundary condition on external boundaries; 
    // reflecting on inherently reflecting boundaries (lox (radial axis),
    // loy, hiy)
    for (int cell = 1; cell <= ncells; cell++)
    {
	layout(cell, 1, 1) = cell;
	layout(cell, 2, 1) = 0; 
	layout(cell, 3, 1) = cell;
	layout(cell, 4, 1) = cell;
	layout(cell, 5, 1) = 0;
	layout(cell, 6, 1) = 0;
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
    if (!soft_equiv(mesh->get_low_x(1),  0.0))  ITFAILS;
    if (!soft_equiv(mesh->get_high_x(1), 1.0))  ITFAILS;
    if (!soft_equiv(mesh->get_low_z(1),  0.0))  ITFAILS;
    if (!soft_equiv(mesh->get_high_z(1), 1.0))  ITFAILS;

    // check that cell midpoints are correct
    if (!soft_equiv(mesh->get_x_midpoint(1), 0.5)) ITFAILS;
    if (!soft_equiv(mesh->get_y_midpoint(1), 0.0)) ITFAILS;
    if (!soft_equiv(mesh->get_z_midpoint(1), 0.5)) ITFAILS;

    // check that the cell dimensions are correct
    if (!soft_equiv(mesh->dim(1, 1), 1.0)) ITFAILS;
    if (!soft_equiv(mesh->dim(2, 1), 1.0)) ITFAILS;
    if (!soft_equiv(mesh->dim(3, 1), 1.0)) ITFAILS;

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
    // on a vacuum boundary
    vector<int> sur_cell(1,1);
    if ( mesh->check_defined_surcells("lox",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
    if ( mesh->check_defined_surcells("loy",sur_cell)) ITFAILS;
    if ( mesh->check_defined_surcells("hiy",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("loz",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("hiz",sur_cell)) ITFAILS;

    // check that the mesh returns the proper cell list (only 1 cell in mesh)
    vector<int> surface_cell_list;
    surface_cell_list = mesh->get_surcells("hir");
    if (surface_cell_list.size() != 1)             ITFAILS;
    if (surface_cell_list[0] != 1)                 ITFAILS;
    surface_cell_list = mesh->get_surcells("loz");
    if (surface_cell_list.size() != 1)             ITFAILS;
    if (surface_cell_list[0] != 1)                 ITFAILS;
    surface_cell_list = mesh->get_surcells("hiz");
    if (surface_cell_list.size() != 1)             ITFAILS;
    if (surface_cell_list[0] != 1)                 ITFAILS;

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

    // >>>> test the mesh's ability to sample positions. <<<<

    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    if (seed != 1234567) ITFAILS;

    // The explicit-value tests has been hand-checked in a weak sense and,
    // therefore, are more of a regression test.
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

    // sample with slope -- explicit-value test: 
    vector<double> slope(3);
    slope[0] = -1.0;
    slope[1] =  0.0;
    slope[2] =  0.0;

    position_sampled = mesh->sample_pos(1, ran_object, slope, 1.0);
    if (position_sampled.size() != 3)                    ITFAILS;
    if (!soft_equiv(position_sampled[0],  0.2633,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[1],-.002651,0.001)) ITFAILS;
    if (!soft_equiv(position_sampled[2],  0.0627,0.001)) ITFAILS;

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

    // with slope=(-1,0,0) and a cell value T^4 of 0.5, x- and z-dimensions
    // should be sampled uniformly
    int num_bins = 100;
    vector<double> xdist(num_bins, 0.0);
    vector<double> zdist(num_bins, 0.0);
    double xavg = 0.0;
    double yavg = 0.0;
    double zavg = 0.0;
    double fabs_y_avg = 0.0;
    int xbin, zbin;
    int num_particles = 10000;

    slope[0] = -1.0;
    slope[1] =  0.0;
    slope[2] =  0.0;
    
    for (int num_p = 0; num_p < num_particles; num_p++)
    {
	position_sampled = mesh->sample_pos(1, ran_object, slope, 0.5);
	
	xbin = position_sampled[0]*num_bins;
	zbin = position_sampled[2]*num_bins;

	xdist[xbin]++;
	zdist[zbin]++;

	xavg += position_sampled[0];
	yavg += position_sampled[1];
	zavg += position_sampled[2];

	fabs_y_avg += std::fabs(position_sampled[1]);
    }

    xavg = xavg/num_particles;
    yavg = yavg/num_particles;
    zavg = zavg/num_particles;
    fabs_y_avg = fabs_y_avg/num_particles;

    double est_std_of_mean = 1.0/std::pow(num_particles, 0.5);

    if (!soft_equiv(xavg, 0.5, 4.0*est_std_of_mean))  ITFAILS;
    if (!soft_equiv(yavg, 0.0, 4.0*est_std_of_mean))  ITFAILS;
    if (!soft_equiv(zavg, 0.5, 4.0*est_std_of_mean))  ITFAILS;

    // what is the average |y| given that x is uniform in (0,1)?  noting
    // that, for any x, half the y's will lie under x/2, thus half the y's
    // are under the line y=x/2, and the other half of the y's are between
    // y=x/2 and y=x.  Integrating, then |y|_avg = int_0^1 (x/2) = 1/4
    if (!soft_equiv(fabs_y_avg, 0.25, 4.0*est_std_of_mean))  ITFAILS;
    
    for (int bin = 0; bin < num_bins; bin++)
    {
	xdist[bin] = xdist[bin] / num_particles * num_bins;
	zdist[bin] = zdist[bin] / num_particles * num_bins;
    }

    int num_xbin_within_1std = 0;
    int num_zbin_within_1std = 0;
    int num_xbin_within_2std = 0;
    int num_zbin_within_2std = 0;
    int num_xbin_within_3std = 0;
    int num_zbin_within_3std = 0;

    double estimated_std = 
	std::pow(static_cast<double>(num_bins)/num_particles,0.5);

    for (int bin = 0; bin < num_bins; bin++)
    {
	if (std::fabs(xdist[bin]-1.0) < estimated_std)
	    num_xbin_within_1std++;
	if (std::fabs(zdist[bin]-1.0) < estimated_std)
	    num_zbin_within_1std++;

	if (std::fabs(xdist[bin]-1.0) < 2.0*estimated_std)
	    num_xbin_within_2std++;
	if (std::fabs(zdist[bin]-1.0) < 2.0*estimated_std)
	    num_zbin_within_2std++;

	if (std::fabs(xdist[bin]-1.0) < 3.0*estimated_std)
	    num_xbin_within_3std++;
	if (std::fabs(zdist[bin]-1.0) < 3.0*estimated_std)
	    num_zbin_within_3std++;
    }
    double coverage_x_1std =
	static_cast<double>(num_xbin_within_1std)/num_bins; 
    double coverage_z_1std =
	static_cast<double>(num_zbin_within_1std)/num_bins; 
    double coverage_x_2std =
	static_cast<double>(num_xbin_within_2std)/num_bins; 
    double coverage_z_2std =
	static_cast<double>(num_zbin_within_2std)/num_bins; 
    double coverage_x_3std =
	static_cast<double>(num_xbin_within_3std)/num_bins; 
    double coverage_z_3std =
	static_cast<double>(num_zbin_within_3std)/num_bins; 

    if ((coverage_x_1std < 0.67) && 
	!soft_equiv(coverage_x_1std, 0.67, 0.1))	ITFAILS; 
    if ((coverage_z_1std < 0.67) && 
	!soft_equiv(coverage_z_1std, 0.67, 0.1))	ITFAILS; 
    if ((coverage_x_2std < 0.95) && 
	!soft_equiv(coverage_x_2std, 0.95, 0.1))	ITFAILS;  
    if ((coverage_z_2std < 0.95) && 
	!soft_equiv(coverage_z_2std, 0.95, 0.1))	ITFAILS;  
    if ((coverage_x_3std < 0.99) && 
	!soft_equiv(coverage_x_3std, 0.99, 0.1))	ITFAILS;  
    if ((coverage_z_3std < 0.99) && 
	!soft_equiv(coverage_z_3std, 0.99, 0.1))	ITFAILS;  

    // <<<<< uniform sampling in volume >>>>>>
    // now, let's sample again, except totally uniform in volume, and just
    // check the averages
    xavg       = 0.0;
    yavg       = 0.0;
    zavg       = 0.0;
    fabs_y_avg = 0.0;

    slope[0] = 0.0;
    slope[1] = 0.0;
    slope[2] = 0.0;
    
    for (int num_p = 0; num_p < num_particles; num_p++)
    {
	position_sampled = mesh->sample_pos(1, ran_object, slope, 0.5);
	
	xavg += position_sampled[0];
	yavg += position_sampled[1];
	zavg += position_sampled[2];

	fabs_y_avg += std::fabs(position_sampled[1]);
    }

    xavg       = xavg/num_particles;
    yavg       = yavg/num_particles;
    zavg       = zavg/num_particles;
    fabs_y_avg = fabs_y_avg/num_particles;

    // now for a sampling that is uniform in volume, we see that, for this
    // one-cell problem of 90 degrees, f(x)=2x and f(|y|)=2(1-y).  Thus,
    // x_avg   = int_0^1(xf(x))dx = int_0^1(2x^2)dx    = 2/3
    // |y|_avg = int_0^1(yf(y))dy = int_0^1(y2(1-y))dy = 1/3

    est_std_of_mean = 1.0/std::pow(num_particles, 0.5);

    if (!soft_equiv(xavg, 2./3., 4.0*est_std_of_mean))  ITFAILS;
    if (!soft_equiv(yavg,   0.0, 4.0*est_std_of_mean))  ITFAILS;
    if (!soft_equiv(zavg,   0.5, 4.0*est_std_of_mean))  ITFAILS;

    if (!soft_equiv(fabs_y_avg, 1./3., 4.0*est_std_of_mean))  ITFAILS;
    

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
// test the RZWedge_Mesh via the RZWedge_Builder

void build_an_RZWedge()
{
    // make a builder from parsing the RZWedge input
    SP<Parser_RZ> parser(new Parser_RZ());
    RZWedge_Builder builder(parser);

    // check some of the RZWedge_Builder properties; before build_Mesh, the
    // coordinate system is two-dimensional.  After build_Mesh, it's 3D XYZ. 

    // test cell region data (same as in OS_Mesh)
    {
	vector<int> regions(6,1);
	regions[3] = 2;
	regions[4] = 2;
	regions[5] = 2;

	if (builder.get_num_regions() != 2)   ITFAILS;
	if (builder.get_regions() != regions) ITFAILS;
    }

    // test zone mapping (same as in OS_Mesh)
    {
	vector<vector<int> > zone(2, vector<int>(3));
	zone[0][0] = 1;
	zone[0][1] = 2;
	zone[0][2] = 3;
	zone[1][0] = 4;
	zone[1][1] = 5;
	zone[1][2] = 6;

	if (builder.get_num_zones() != 2)            ITFAILS;
	if (builder.get_cells_in_zone(1) != zone[0]) ITFAILS;
	if (builder.get_cells_in_zone(2) != zone[1]) ITFAILS;
    }

    // build an RZWedge mesh
    SP<RZWedge_Mesh> mesh = builder.build_Mesh();

    // check defined surface source cells
    {
	vector<vector<int> > ss(2);
	ss[0].resize(2);
	ss[1].resize(2);
	ss[0][0] = 3;
	ss[0][1] = 6;
	ss[1][0] = 1;
	ss[1][1] = 2;

	if (builder.get_defined_surcells() != ss) ITFAILS;

	vector<string> ssp(2);
	ssp[0] = "hir";
	ssp[1] = "loz";
	
	if (builder.get_ss_pos() != ssp) ITFAILS;
    }

    // check zone mapper
    {
	vector<int> zone_field(2);
	zone_field[0] = 1000;
	zone_field[1] = 1001;
	vector<int> cell_field = builder.zone_cell_mapper(zone_field);

	if (cell_field.size() != mesh->num_cells()) ITFAILS;

	if (cell_field[0] != 1000) ITFAILS;
	if (cell_field[1] != 1000) ITFAILS;
	if (cell_field[2] != 1000) ITFAILS;
	if (cell_field[3] != 1001) ITFAILS;
	if (cell_field[4] != 1001) ITFAILS;
	if (cell_field[5] != 1001) ITFAILS;
    }

    //     <<<< now do some checks on the mesh itself >>>>>
    // (this is a different mesh than the simple_one_cell_RZWedge)

    // check that the mesh returns proper coordinate system information
    if (mesh->get_Coord().get_dim() != 3)    ITFAILS;
    if (mesh->get_SPCoord()->get_dim() != 3) ITFAILS;
    if (mesh->get_Coord().get_Coord()    != std::string("xyz")) ITFAILS;
    if (mesh->get_SPCoord()->get_Coord() != std::string("xyz")) ITFAILS;

    // check that mesh returns the proper number of cells
    if (mesh->num_cells() != 6) ITFAILS;

    // calculate sqrt(pi)/2 -- it will come in handy!
    double sqrtpi_2 = 0.5 * std::pow(pi,0.5);

    // check that x- and z-extents are correct for the 1st and last cells
    if (!soft_equiv(mesh->get_low_x(1), 0.0))       ITFAILS;
    if (!soft_equiv(mesh->get_high_x(1), sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->get_low_z(1), -1.0))      ITFAILS;
    if (!soft_equiv(mesh->get_high_z(1), 1.0))      ITFAILS;

    if (!soft_equiv(mesh->get_low_x(6), 2.0*sqrtpi_2))  ITFAILS;
    if (!soft_equiv(mesh->get_high_x(6), 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->get_low_z(6), 1.0))           ITFAILS;
    if (!soft_equiv(mesh->get_high_z(6), 3.0))          ITFAILS;

    // check that cell midpoints are correct for 1st and last cells
    if (!soft_equiv(mesh->get_x_midpoint(1), 0.5*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->get_y_midpoint(1), 0.0))          ITFAILS;
    if (!soft_equiv(mesh->get_z_midpoint(1), 0.0))          ITFAILS;

    if (!soft_equiv(mesh->get_x_midpoint(6), 2.5*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->get_y_midpoint(6), 0.0))          ITFAILS;
    if (!soft_equiv(mesh->get_z_midpoint(6), 2.0))          ITFAILS;

    // check that the cell dimensions are correct for 1st and last cells
    if (!soft_equiv(mesh->dim(1,1), sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->dim(2,1), sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->dim(3,1), 2.0))      ITFAILS;

    if (!soft_equiv(mesh->dim(1,6), sqrtpi_2))     ITFAILS;
    if (!soft_equiv(mesh->dim(2,6), 5.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(mesh->dim(3,6), 2.0))          ITFAILS;

    // check that graphics cell type is correct
    vector<int> cell_type = mesh->get_cell_types();
    if (cell_type.size() != 6)                               ITFAILS;
    for (int i = 0; i < 6; i++) 
	if (cell_type[0] != rtt_viz::eight_node_hexahedron ) ITFAILS;

    // check that mesh returns the proper point coordinates (1st & last cell)
    vector<vector<double> > point_coords = mesh->get_point_coord();
    if (point_coords.size() != 48)       ITFAILS;
    for (int p = 0; p < 48; p++) 
	if (point_coords[p].size() != 3) ITFAILS;

    //   first cell, low z face, clockwise from the lower left, facing +z
    if (!soft_equiv(point_coords[0][0], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[0][1], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[0][2],-1.0))      ITFAILS;
    if (!soft_equiv(point_coords[1][0], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[1][1],-sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[1][2], -1.0))     ITFAILS;
    if (!soft_equiv(point_coords[2][0], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[2][1], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[2][2], -1.0))     ITFAILS;
    if (!soft_equiv(point_coords[3][0], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[3][1], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[3][2], -1.0))     ITFAILS;

    //   first cell, high z face, clockwise from the lower left, facing +z 
    if (!soft_equiv(point_coords[4][0], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[4][1], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[4][2], 1.0))      ITFAILS;
    if (!soft_equiv(point_coords[5][0], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[5][1],-sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[5][2], 1.0))      ITFAILS;
    if (!soft_equiv(point_coords[6][0], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[6][1], sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[6][2], 1.0))      ITFAILS;
    if (!soft_equiv(point_coords[7][0], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[7][1], 0.0))      ITFAILS;
    if (!soft_equiv(point_coords[7][2], 1.0))      ITFAILS;

    //   last cell, low z face, clockwise from the lower left, facing +z
    if (!soft_equiv(point_coords[40][0], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[40][1],-2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[40][2], 1.0))          ITFAILS;
    if (!soft_equiv(point_coords[41][0], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[41][1],-3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[41][2], 1.0))          ITFAILS;
    if (!soft_equiv(point_coords[42][0], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[42][1], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[42][2], 1.0))          ITFAILS;
    if (!soft_equiv(point_coords[43][0], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[43][1], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[43][2], 1.0))          ITFAILS;

    //   last cell, high z face, clockwise from the lower left, facing +z 
    if (!soft_equiv(point_coords[44][0], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[44][1],-2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[44][2], 3.0))          ITFAILS;
    if (!soft_equiv(point_coords[45][0], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[45][1],-3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[45][2], 3.0))          ITFAILS;
    if (!soft_equiv(point_coords[46][0], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[46][1], 3.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[46][2], 3.0))          ITFAILS;
    if (!soft_equiv(point_coords[47][0], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[47][1], 2.0*sqrtpi_2)) ITFAILS;
    if (!soft_equiv(point_coords[47][2], 3.0))          ITFAILS;

    // check that mesh returns the correct "next cell"
    // wedge faces (low and high y) should be reflecting
    for (int c = 1; c <= 6; c++)
    {
	if (mesh->next_cell(c, 3) != c) ITFAILS;
	if (mesh->next_cell(c, 4) != c) ITFAILS;
    }

    if (mesh->next_cell(1, 1) != 1) ITFAILS;
    if (mesh->next_cell(1, 2) != 2) ITFAILS;
    if (mesh->next_cell(1, 5) != 0) ITFAILS;
    if (mesh->next_cell(1, 6) != 4) ITFAILS;
    if (mesh->next_cell(2, 1) != 1) ITFAILS;
    if (mesh->next_cell(2, 2) != 3) ITFAILS;
    if (mesh->next_cell(2, 5) != 0) ITFAILS;
    if (mesh->next_cell(2, 6) != 5) ITFAILS;
    if (mesh->next_cell(3, 1) != 2) ITFAILS;
    if (mesh->next_cell(3, 2) != 0) ITFAILS;
    if (mesh->next_cell(3, 5) != 0) ITFAILS;
    if (mesh->next_cell(3, 6) != 6) ITFAILS;
    if (mesh->next_cell(4, 1) != 4) ITFAILS;
    if (mesh->next_cell(4, 2) != 5) ITFAILS;
    if (mesh->next_cell(4, 5) != 1) ITFAILS;
    if (mesh->next_cell(4, 6) != 0) ITFAILS;
    if (mesh->next_cell(5, 1) != 4) ITFAILS;
    if (mesh->next_cell(5, 2) != 6) ITFAILS;
    if (mesh->next_cell(5, 5) != 2) ITFAILS;
    if (mesh->next_cell(5, 6) != 0) ITFAILS;
    if (mesh->next_cell(6, 1) != 5) ITFAILS;
    if (mesh->next_cell(6, 2) != 0) ITFAILS;
    if (mesh->next_cell(6, 5) != 3) ITFAILS;
    if (mesh->next_cell(6, 6) != 0) ITFAILS;

    // check that mesh can locate a cell given a position
    vector<double> position(3, 0.0);
    position[0] = 0.5;
    position[2] = 0.5;
    if (mesh->get_cell(position) != 1)  ITFAILS;

    position[0] = 100.0;
    if (mesh->get_cell(position) != -1) ITFAILS;

    position[0] = 1.5;
    position[1] = 1.0;
    position[2] = 2.0;
    if (mesh->get_cell(position) != 5)  ITFAILS;


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

    // now see that all other cells have the same normals
    for (int c = 1; c <= 6; c++)
	for (int f = 1; f <= 6; f++)
	{
	    normal = mesh->get_normal(1,f);
	    vector<double> check = mesh->get_normal(c,f);
	    for (int d = 0; d < 3; d++)
		if (!soft_equiv(normal[d],check[d]))    ITFAILS;
	}
		    

    vector<double> normal_in(3, 0.0);
    for (int cell = 1; cell <= 6; cell++)
	for (int face = 1; face <= 6; face++)
	{
	    normal = mesh->get_normal(1,face);
	    normal_in = mesh->get_normal_in(1,face);
	    for (int dim = 0; dim < 3; dim++)
		if (!soft_equiv(normal[dim],-normal_in[dim])) ITFAILS;
	}
    
    // does the mesh return the correct cell volume?
    for (int cell = 1; cell <= 3; cell++)
    {
	if (!soft_equiv(mesh->volume(cell),  (2*cell-1)*pi*0.5)) ITFAILS;
	if (!soft_equiv(mesh->volume(cell+3),(2*cell-1)*pi*0.5)) ITFAILS;
    }


    // does the mesh return the correct face areas?
    if (!soft_equiv(mesh->face_area(1,1),0.0))             ITFAILS;
    if (!soft_equiv(mesh->face_area(1,2),pow(4.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(1,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(1,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(1,5),pi*0.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(1,6),pi*0.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(4,1),0.0))             ITFAILS;
    if (!soft_equiv(mesh->face_area(4,2),pow(4.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(4,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(4,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(4,5),pi*0.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(4,6),pi*0.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(2,1),pow(4.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(2,2),pow(16.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(2,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(2,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(2,5),pi*0.75))         ITFAILS;
    if (!soft_equiv(mesh->face_area(2,6),pi*0.75))         ITFAILS;
    if (!soft_equiv(mesh->face_area(5,1),pow(4.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(5,2),pow(16.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(5,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(5,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(5,5),pi*0.75))         ITFAILS;
    if (!soft_equiv(mesh->face_area(5,6),pi*0.75))         ITFAILS;
    if (!soft_equiv(mesh->face_area(3,1),pow(16.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(3,2),pow(36.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(3,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(3,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(3,5),pi*1.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(3,6),pi*1.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(6,1),pow(16.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(6,2),pow(36.0*pi,0.5)))ITFAILS;
    if (!soft_equiv(mesh->face_area(6,3),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(6,4),pow(2.0*pi,0.5))) ITFAILS;
    if (!soft_equiv(mesh->face_area(6,5),pi*1.25))         ITFAILS;
    if (!soft_equiv(mesh->face_area(6,6),pi*1.25))         ITFAILS;

    // check that the mesh returns correctly sized vertices.
    vector<vector<double> > vertices;
    for (int cell = 0; cell < 6; cell++)
    {
	for (int face = 0; face < 6; face++)
	{
	    vertices = mesh->get_vertices(cell+1,face+1);
	    
	    if (vertices.size() != 3)    ITFAILS;  // 3 dimensions
	    if (vertices[0].size() != 4) ITFAILS;  // 4 nodes per face
	    if (vertices[1].size() != 4) ITFAILS;  // 4 nodes per face
	    if (vertices[2].size() != 4) ITFAILS;  // 4 nodes per face
	}

	vertices = mesh->get_vertices(cell+1);
	if (vertices.size() != 3)    ITFAILS;  // 3 dimensions
	if (vertices[0].size() != 8) ITFAILS;  // 8 nodes per cell
	if (vertices[1].size() != 8) ITFAILS;  // 8 nodes per cell
	if (vertices[2].size() != 8) ITFAILS;  // 8 nodes per cell
    }

    // spot check a few vertices
    vertices = mesh->get_vertices(6,5);
    if (!soft_equiv(vertices[0][0], pow(pi,0.5))) ITFAILS;
    if (!soft_equiv(vertices[1][0],-pow(pi,0.5))) ITFAILS;
    if (!soft_equiv(vertices[2][0], 1.0))         ITFAILS;
    vertices = mesh->get_vertices(2,2);
    if (!soft_equiv(vertices[0][2], pow(pi,0.5))) ITFAILS;
    if (!soft_equiv(vertices[1][2], pow(pi,0.5))) ITFAILS;
    if (!soft_equiv(vertices[2][2], 1.0))         ITFAILS;
    

    // check that the mesh returns the correct face number
    if (mesh->get_bndface("lox", 1) != 1) ITFAILS;
    if (mesh->get_bndface("hix", 1) != 2) ITFAILS;
    if (mesh->get_bndface("lor", 1) != 1) ITFAILS;
    if (mesh->get_bndface("hir", 1) != 2) ITFAILS;
    if (mesh->get_bndface("loy", 1) != 3) ITFAILS;
    if (mesh->get_bndface("hiy", 1) != 4) ITFAILS;
    if (mesh->get_bndface("loz", 1) != 5) ITFAILS;
    if (mesh->get_bndface("hiz", 1) != 6) ITFAILS;

    // check that surface source cells would actually be on a vacuum boundary
    vector<int> sur_cell;
    sur_cell.resize(2);
    sur_cell[0] = 3;
    sur_cell[1] = 6;
    if (!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
    if (!mesh->check_defined_surcells("hir",sur_cell)) ITFAILS;
    sur_cell[0] = 1;
    sur_cell[1] = 4;
    if (mesh->check_defined_surcells("lox",sur_cell))  ITFAILS;
    if (mesh->check_defined_surcells("lor",sur_cell))  ITFAILS;
    sur_cell.resize(3);
    sur_cell[0] = 1;
    sur_cell[1] = 2;
    sur_cell[2] = 3;
    if (!mesh->check_defined_surcells("loz",sur_cell)) ITFAILS;
    sur_cell[0] = 4;
    sur_cell[1] = 5;
    sur_cell[2] = 6;
    if (!mesh->check_defined_surcells("hiz",sur_cell)) ITFAILS;
    sur_cell.resize(6);
    for (int c = 0; c < 6; c++)
	sur_cell[c] = c+1;
    if ( mesh->check_defined_surcells("loy",sur_cell)) ITFAILS;
    if ( mesh->check_defined_surcells("hiy",sur_cell)) ITFAILS;

    // check that the mesh returns the proper cell list 
    vector<int> surface_cell_list;
    surface_cell_list = mesh->get_surcells("hir");
    if (surface_cell_list.size() != 2)             ITFAILS;
    if (surface_cell_list[0] != 3)                 ITFAILS;
    if (surface_cell_list[1] != 6)                 ITFAILS;
    surface_cell_list = mesh->get_surcells("loz");
    if (surface_cell_list.size() != 3)             ITFAILS;
    if (surface_cell_list[0] != 1)                 ITFAILS;
    if (surface_cell_list[1] != 2)                 ITFAILS;
    if (surface_cell_list[2] != 3)                 ITFAILS;
    surface_cell_list = mesh->get_surcells("hiz");
    if (surface_cell_list.size() != 3)             ITFAILS;
    if (surface_cell_list[0] != 4)                 ITFAILS;
    if (surface_cell_list[1] != 5)                 ITFAILS;
    if (surface_cell_list[2] != 6)                 ITFAILS;

    // check that the neighbors are correct
    vector<int> neighbors;
    neighbors = mesh->get_neighbors(1);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 1)     ITFAILS;
    if (neighbors[1] != 2)     ITFAILS;
    if (neighbors[2] != 1)     ITFAILS;
    if (neighbors[3] != 1)     ITFAILS;
    if (neighbors[4] != 0)     ITFAILS;
    if (neighbors[5] != 4)     ITFAILS;
    neighbors = mesh->get_neighbors(2);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 1)     ITFAILS;
    if (neighbors[1] != 3)     ITFAILS;
    if (neighbors[2] != 2)     ITFAILS;
    if (neighbors[3] != 2)     ITFAILS;
    if (neighbors[4] != 0)     ITFAILS;
    if (neighbors[5] != 5)     ITFAILS;
    neighbors = mesh->get_neighbors(3);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 2)     ITFAILS;
    if (neighbors[1] != 0)     ITFAILS;
    if (neighbors[2] != 3)     ITFAILS;
    if (neighbors[3] != 3)     ITFAILS;
    if (neighbors[4] != 0)     ITFAILS;
    if (neighbors[5] != 6)     ITFAILS;
    neighbors = mesh->get_neighbors(4);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 4)     ITFAILS;
    if (neighbors[1] != 5)     ITFAILS;
    if (neighbors[2] != 4)     ITFAILS;
    if (neighbors[3] != 4)     ITFAILS;
    if (neighbors[4] != 1)     ITFAILS;
    if (neighbors[5] != 0)     ITFAILS;
    neighbors = mesh->get_neighbors(5);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 4)     ITFAILS;
    if (neighbors[1] != 6)     ITFAILS;
    if (neighbors[2] != 5)     ITFAILS;
    if (neighbors[3] != 5)     ITFAILS;
    if (neighbors[4] != 2)     ITFAILS;
    if (neighbors[5] != 0)     ITFAILS;
    neighbors = mesh->get_neighbors(6);
    if (neighbors.size() != 6) ITFAILS;
    if (neighbors[0] != 5)     ITFAILS;
    if (neighbors[1] != 0)     ITFAILS;
    if (neighbors[2] != 6)     ITFAILS;
    if (neighbors[3] != 6)     ITFAILS;
    if (neighbors[4] != 3)     ITFAILS;
    if (neighbors[5] != 0)     ITFAILS;

    // make sure the mesh thinks it's not a submesh
    if (!mesh->full_Mesh()) ITFAILS;

}

//---------------------------------------------------------------------------//
// Test AMR features of the RZWedge_Mesh

void test_AMR()
{
    // make a RZWedge_AMR mesh
    SP<RZWedge_Mesh> mesh = make_RZWedge_Mesh_AMR(10);
    if (mesh->num_cells() != 12) ITFAILS;

    // check neighbors function
    vector<int> nbors;

    // cell 1
    {
	nbors.resize(7);
	nbors[0] = 1;
	nbors[1] = 7;
	nbors[2] = 1;
	nbors[3] = 1;
	nbors[4] = 1;
	nbors[5] = 2;
	nbors[6] = 4;

	if (mesh->get_neighbors(1) != nbors) ITFAILS;
    }

    // cell 3
    {
	nbors.resize(6);
	nbors[0] = 3;
	nbors[1] = 5;
	nbors[2] = 3;
	nbors[3] = 3;
	nbors[4] = 2;
	nbors[5] = 6;

	if (mesh->get_neighbors(3) != nbors) ITFAILS;
    }

    // cell 6
    {
	nbors.resize(8);
	nbors[0] = 6;
	nbors[1] = 9;
	nbors[2] = 10;
	nbors[3] = 6;
	nbors[4] = 6;
	nbors[5] = 3;
	nbors[6] = 5;
	nbors[7] = 0;

	if (mesh->get_neighbors(6) != nbors) ITFAILS;
    }

    // cell 7
    {
	nbors.resize(6);
	nbors[0] = 1;
	nbors[1] = 0;
	nbors[2] = 7;
	nbors[3] = 7;
	nbors[4] = 7;
	nbors[5] = 8;

	if (mesh->get_neighbors(7) != nbors) ITFAILS;
    }

    // cell 8
    {
	nbors.resize(8);
	nbors[0] = 4;
	nbors[1] = 5;
	nbors[2] = 0;
	nbors[3] = 8;
	nbors[4] = 8;
	nbors[5] = 7;
	nbors[6] = 9;
	nbors[7] = 11;

	if (mesh->get_neighbors(8) != nbors) ITFAILS;
    }

    // cell 11
    {
	nbors.resize(6);
	nbors[0] = 9;
	nbors[1] = 0;
	nbors[2] = 11;
	nbors[3] = 11;
	nbors[4] = 8;
	nbors[5] = 12;

	if (mesh->get_neighbors(11) != nbors) ITFAILS;
    }

    // check next_cell function
    
    // position vector (function calls without r checked above)
    vector<double> r(3);

    // cell 1
    {
	if (mesh->next_cell(1, 1, r) != 1) ITFAILS;
	if (mesh->next_cell(1, 2, r) != 7) ITFAILS;
	if (mesh->next_cell(1, 3, r) != 1) ITFAILS;
	if (mesh->next_cell(1, 4, r) != 1) ITFAILS;
	if (mesh->next_cell(1, 5, r) != 1) ITFAILS;

	r[0] = 0.501;
	r[2] = 1.0;
	if (mesh->next_cell(1, 6, r) != 4) ITFAILS;

	r[0] = 0.499;
	if (mesh->next_cell(1, 6, r) != 2) ITFAILS;

	r[0] = 0.500;
	if (mesh->next_cell(1, 6, r) != 2) ITFAILS;

	r[0] = 1.0000001;  // valid to 6 digits 
	if (mesh->next_cell(1, 6, r) != 4) ITFAILS;
	
	r[0] = -0.0000001; // valid to 6 digits
	if (mesh->next_cell(1, 6, r) != 2) ITFAILS;
    }

    // cell 3
    {
	// only 1 cell across faces; r is not used
	if (mesh->next_cell(3, 1, r) != 3) ITFAILS;
	if (mesh->next_cell(3, 2, r) != 5) ITFAILS;
	if (mesh->next_cell(3, 3, r) != 3) ITFAILS;
	if (mesh->next_cell(3, 4, r) != 3) ITFAILS;
	if (mesh->next_cell(3, 5, r) != 2) ITFAILS;
	if (mesh->next_cell(3, 6, r) != 6) ITFAILS;
    }

    // cell 6
    {
	if (mesh->next_cell(6, 1, r) != 6)  ITFAILS;

	r[0] = 1.0000001; // must be on surface within 1.e-6
	r[2] = 2.499;
	if (mesh->next_cell(6, 2, r) != 9)  ITFAILS;

	r[2] = 2.501;
	if (mesh->next_cell(6, 2, r) != 10) ITFAILS;

	r[2] = 2.500;
	if (mesh->next_cell(6, 2, r) != 9)  ITFAILS;

	r[2] = 3.000001;
	if (mesh->next_cell(6, 2, r) != 10) ITFAILS;

	if (mesh->next_cell(6, 3, r) != 6)  ITFAILS;
	if (mesh->next_cell(6, 4, r) != 6)  ITFAILS;

	r[0] = 0.499;
	r[2] = 2.0 - .0000001;
	if (mesh->next_cell(6, 5, r) != 3)  ITFAILS;
	
	r[0] = -0.0000001;
	if (mesh->next_cell(6, 5, r) != 3)  ITFAILS;
	
	r[0] = .50001;
	if (mesh->next_cell(6, 5, r) != 5)  ITFAILS;

	r[0] = .5000;
	if (mesh->next_cell(6, 5, r) != 3)  ITFAILS;
	
	if (mesh->next_cell(6, 6, r) != 0)  ITFAILS;
    }

    // cell 8
    {
	r[0] = 1.0 - 0.0000001; // must be on surface within 1.e-6
	r[2] = 1.499;
	if (mesh->next_cell(8, 1, r) != 4)  ITFAILS;

	r[2] = 1.501;
	if (mesh->next_cell(8, 1, r) != 5)  ITFAILS;

	r[2] = 2.0000001;
	if (mesh->next_cell(8, 1, r) != 5)  ITFAILS;

	if (mesh->next_cell(8, 2, r) != 0)  ITFAILS;
	if (mesh->next_cell(8, 3, r) != 8)  ITFAILS;
	if (mesh->next_cell(8, 4, r) != 8)  ITFAILS;
	if (mesh->next_cell(8, 5, r) != 7)  ITFAILS;

	r[0] = 1.499;
	r[2] = 2.0000001;
	if (mesh->next_cell(8, 6, r) != 9)  ITFAILS;
	
	r[0] = 1.0 - 0.0000001;
	if (mesh->next_cell(8, 6, r) != 9)  ITFAILS;
	
	r[0] = 1.50001;
	if (mesh->next_cell(8, 6, r) != 11) ITFAILS;

	r[0] = 1.5000;
	if (mesh->next_cell(8, 6, r) != 9)  ITFAILS;
    }

    // check surface cell functions
    vector<int> sc;

    // high r
    {
	sc.resize(4);
	sc[0] = 7;
	sc[1] = 8;
	sc[2] = 11;
	sc[3] = 12;

	if (mesh->get_surcells("hix") != sc) ITFAILS;
    }

    // low z
    {
	sc.resize(2);
	sc[0] = 1;
	sc[1] = 7;

	if (mesh->get_surcells("loz") != sc) ITFAILS;
    }

    // high z
    {
	sc.resize(3);
	sc[0] = 6;
	sc[1] = 10;
	sc[2] = 12;

	if (mesh->get_surcells("hiz") != sc) ITFAILS;
    }

    // check in cell function
    
    // cell 5
    {
	r[0] = .5 - 1.e-13;
	r[1] = 0.0;
	r[2] = 1.5;

	if (!mesh->in_cell(5, r)) ITFAILS;
    }

    // cell 8
    {
	r[0] = 2.0;
	r[1] = 0.0;
	r[2] = 2.0 + 1.e-13;
	
	if (!mesh->in_cell(8, r)) ITFAILS;
    }
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

	// test the RZWedge_Mesh Builder
	build_an_RZWedge();

	// test AMR features of the mesh
	test_AMR();
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
