//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSphyramid_Mesh.cc
 * \author Jeffery Densmore (again, stolen from tstRZWedge_Mesh.cc)
 * \date   Mon Nov 10 9:51:00 2003
 * \brief  
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../Sphyramid_Mesh.hh"
#include "../XYZCoord_sys.hh"
#include "../Layout.hh"
#include "../Sphyramid_Builder.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "viz/Ensight_Translator.hh"
#include "rng/Rnd_Control.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void simple_one_cell_Sphyramid()
{
    using rtt_dsxx::SP;
    using rtt_mc::XYZCoord_sys;
    using rtt_mc::Layout;
    using rtt_mc::Sphyramid_Mesh;
    using std::string;
    using std::vector;
    using rtt_dsxx::soft_equiv;
    using std::tan;
    using std::sin;
    using std::cos;
    using rtt_rng::Rnd_Control;
    using rtt_rng::Sprng;
    using std::fabs;
    using std::pair;

    // >>> setup one cell mesh <<<

    // one-cell problem
    int ncells = 1 ;
    
    // build an XYZ coordinate system class
    SP<XYZCoord_sys> coord(new XYZCoord_sys());

    // build a layout of size ncells 
    Layout layout;
    layout.set_size(ncells);

    // set number of faces per cell (6)
    for (int cell = 1; cell <= ncells; cell++)
    {
	layout.set_size(cell, 6);
    }
    // set vacuum boundary condition on external boundaries;
    // reflecting on inherently reflecting boundaries
    for (int cell = 1; cell <= ncells; cell++)
    {
	layout(cell,1)=cell;
	layout(cell,2)=0;
	layout(cell,3)=cell;
	layout(cell,4)=cell;
	layout(cell,5)=cell;
	layout(cell,6)=cell;
    }

    // set x cell extents
    vector< vector<double> > cell_x_extents;
    cell_x_extents.resize(ncells);
    for (int cell = 1;cell <= ncells; cell++)
    {
	cell_x_extents[cell-1].resize(2);
    }
    for (int cell = 1; cell <= ncells; cell++)
    {
	cell_x_extents[cell-1][0]=0.0; // low x
	cell_x_extents[cell-1][1]=1.0; // high x
    }

    // set Sphyramid angle
    double beta_radians=0.725153345774; // alpha_degrees=45
    
    // calculate angle data, you'll need it
    double tan_beta = tan(beta_radians);
    double cos_beta = cos(beta_radians);
    double sin_beta = sin(beta_radians);

    // build a mesh object
    SP<Sphyramid_Mesh> mesh(new Sphyramid_Mesh(coord, layout, cell_x_extents, 
					   beta_radians));
    
    // check that the mesh returns proper coordinate system information
    if (mesh->get_Coord().get_dim()      != 3)             ITFAILS;
    if (mesh->get_SPCoord()->get_dim()   != 3)             ITFAILS;
    if (mesh->get_Coord().get_Coord()    != string("xyz")) ITFAILS;
    if (mesh->get_SPCoord()->get_Coord() != string("xyz")) ITFAILS;

    // check that mesh returns the proper number of cells
    if (mesh->num_cells() != ncells)  ITFAILS;

    // check that x extents are correct
    if (!soft_equiv(mesh->get_low_x(1),  0.0))  ITFAILS;
    if (!soft_equiv(mesh->get_high_x(1), 1.0))  ITFAILS;

    // check that cell midpoints are correct
    if (!soft_equiv(mesh->get_x_midpoint(1), 0.5)) ITFAILS;
    if (!soft_equiv(mesh->get_y_midpoint(1), 0.0)) ITFAILS;
    if (!soft_equiv(mesh->get_z_midpoint(1), 0.0)) ITFAILS;

    // check that the cell dimensions are correct
    if (!soft_equiv(mesh->dim(1,1), 1.0))       ITFAILS;
    if (!soft_equiv(mesh->dim(2,1), tan_beta))  ITFAILS;
    if (!soft_equiv(mesh->dim(3,1), tan_beta))  ITFAILS;

    // check that graphics cell type is correct
    {
	vector<int> cell_type = mesh->get_cell_types();
	if (cell_type.size() != 1)                          ITFAILS;
	if (cell_type[0] != rtt_viz::eight_node_hexahedron) ITFAILS;
    }

    // check that the mesh returns the proper point coordinates
    {
	vector< vector<double> > point_coords = mesh->get_point_coord();
	if (point_coords.size()    != 8)               ITFAILS;
	if (point_coords[0].size() != 3)               ITFAILS;
	if (point_coords[1].size() != 3)               ITFAILS;
	if (point_coords[2].size() != 3)               ITFAILS;
	if (point_coords[3].size() != 3)               ITFAILS;
	if (point_coords[4].size() != 3)               ITFAILS;
	if (point_coords[5].size() != 3)               ITFAILS;
	if (point_coords[6].size() != 3)               ITFAILS;
	if (point_coords[7].size() != 3)               ITFAILS;

	if (!soft_equiv(point_coords[0][0], 0.0))      ITFAILS;
	if (!soft_equiv(point_coords[0][1], 0.0))      ITFAILS;
	if (!soft_equiv(point_coords[0][2], 0.0))      ITFAILS;

        if (!soft_equiv(point_coords[1][0], 1.0))      ITFAILS;
        if (!soft_equiv(point_coords[1][1],-tan_beta)) ITFAILS;
        if (!soft_equiv(point_coords[1][2],-tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[2][0], 1.0))      ITFAILS;
        if (!soft_equiv(point_coords[2][1], tan_beta)) ITFAILS;
        if (!soft_equiv(point_coords[2][2],-tan_beta)) ITFAILS;

        if (!soft_equiv(point_coords[3][0], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[3][1], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[3][2], 0.0))      ITFAILS;

        if (!soft_equiv(point_coords[4][0], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[4][1], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[4][2], 0.0))      ITFAILS;

        if (!soft_equiv(point_coords[5][0], 1.0))      ITFAILS;
        if (!soft_equiv(point_coords[5][1],-tan_beta)) ITFAILS;
        if (!soft_equiv(point_coords[5][2], tan_beta)) ITFAILS;

        if (!soft_equiv(point_coords[6][0], 1.0))      ITFAILS;
        if (!soft_equiv(point_coords[6][1], tan_beta)) ITFAILS;
        if (!soft_equiv(point_coords[6][2], tan_beta)) ITFAILS;

        if (!soft_equiv(point_coords[7][0], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[7][1], 0.0))      ITFAILS;
        if (!soft_equiv(point_coords[7][2], 0.0))      ITFAILS;
    }

    // check that the mesh returns the correct "next cell"
    if (mesh->next_cell(1,1) != 1) ITFAILS;
    if (mesh->next_cell(1,2) != 0) ITFAILS;
    if (mesh->next_cell(1,3) != 1) ITFAILS;
    if (mesh->next_cell(1,4) != 1) ITFAILS;
    if (mesh->next_cell(1,5) != 1) ITFAILS;
    if (mesh->next_cell(1,6) != 1) ITFAILS;

    // check that mesh can locate a cell give a position
    {
	vector<double> position(3,0.0);
	position[0]=0.5;
	if (mesh->get_cell(position) != 1) ITFAILS;

	position[0]=100;
	if (mesh->get_cell(position) !=-1) ITFAILS;
    }

    // check that the face normals are correct
    {
	vector<double> normal(3,0.0);

	normal = mesh->get_normal(1,1);
	if (!soft_equiv(normal[0],-1.0))      ITFAILS;
	if (!soft_equiv(normal[1], 0.0))      ITFAILS;
	if (!soft_equiv(normal[2], 0.0))      ITFAILS;
    
	normal = mesh->get_normal(1,2);
	if (!soft_equiv(normal[0], 1.0))      ITFAILS;
	if (!soft_equiv(normal[1], 0.0))      ITFAILS;
	if (!soft_equiv(normal[2], 0.0))      ITFAILS;

	normal = mesh->get_normal(1,3);
	if (!soft_equiv(normal[0],-sin_beta)) ITFAILS;
	if (!soft_equiv(normal[1],-cos_beta)) ITFAILS;
	if (!soft_equiv(normal[2], 0.0))      ITFAILS;
    
	normal = mesh->get_normal(1,4);
	if (!soft_equiv(normal[0],-sin_beta)) ITFAILS;
	if (!soft_equiv(normal[1], cos_beta)) ITFAILS;
	if (!soft_equiv(normal[2], 0.0))      ITFAILS;
    
	normal = mesh->get_normal(1,5);
	if (!soft_equiv(normal[0],-sin_beta)) ITFAILS;
        if (!soft_equiv(normal[1], 0.0))      ITFAILS;
	if (!soft_equiv(normal[2],-cos_beta)) ITFAILS;
    
	normal=mesh->get_normal(1,6);
	if (!soft_equiv(normal[0],-sin_beta)) ITFAILS;
	if (!soft_equiv(normal[1], 0.0))      ITFAILS;
	if (!soft_equiv(normal[2], cos_beta)) ITFAILS;

	vector<double> normal_in(3, 0.0);
    
	for (int face = 1; face <= 6; face++)
	{
	    normal    = mesh->get_normal(1,face);
	    normal_in = mesh->get_normal_in(1,face);
	    for (int dim = 0; dim < 3; dim++)
	    {
		if (!soft_equiv(normal[dim],-normal_in[dim])) ITFAILS;
	    }
	}
    }

    // does the mesh return the correct cell volume
    {
	double cell_volume = 4./3.*tan_beta*tan_beta;
	if (!soft_equiv(mesh->volume(1), cell_volume))          ITFAILS;
	if (!soft_equiv(mesh->get_total_volume(), cell_volume)) ITFAILS;
    }



    // does the mesh return the correct face areas
    {
	double area;
	
	area = 0.0;
	if (!soft_equiv(mesh->face_area(1,1), area)) ITFAILS;
    
	area = 4.*tan_beta*tan_beta;
	if (!soft_equiv(mesh->face_area(1,2), area)) ITFAILS;
    
	area = tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(1,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,6), area)) ITFAILS;
    }

    // check that the mesh returns the correct vertices
    {
	// first, set reference vertices
	vector< vector<double> > ref_vert;
	ref_vert.resize(3);
	for(int dim = 0; dim < 3; dim++)
	{
	    ref_vert[dim].resize(8);
	    for(int node = 0; node < 8; node++)
	    {
		ref_vert[dim][node]=0.0;
	    }
	}
	//node 1 (lox, loy, loz) - all zeros
	//node 2 (hix, loy, loz)
	ref_vert[0][1] = 1.0;
	ref_vert[1][1] =-tan_beta;
	ref_vert[2][1] =-tan_beta;
	//node 3 (hix,hiy,loz)
	ref_vert[0][2] = 1.0;
	ref_vert[1][2] = tan_beta;
	ref_vert[2][2] =-tan_beta;
	//node 4 (lox,hiy,loz) - all zeros
	//node 5 (lox, loy, hiz) - all zeros
	//node 6 (hix, loy, hiz)
	ref_vert[0][5] = 1.0;
	ref_vert[1][5] =-tan_beta;
	ref_vert[2][5] = tan_beta;
	//node 7 (hix,hiy,hiz)
	ref_vert[0][6] = 1.0;
	ref_vert[1][6] = tan_beta;
	ref_vert[2][6] = tan_beta;
	//node 8 (lox,hiy,hiz) - all zeros

	vector<vector<int> > node_num(6,vector<int>(0));
	for (int face = 0; face < 6; face++)
	{
	    node_num[face].resize(4);
	}
	//face 1 low x
	node_num[0][0] = 1;
	node_num[0][1] = 4;
	node_num[0][2] = 8;
	node_num[0][3] = 5;
	//face 2 high x
	node_num[1][0] = 2;
	node_num[1][1] = 3;
	node_num[1][2] = 7;
	node_num[1][3] = 6;
	//face 3 low y
	node_num[2][0] = 1;
	node_num[2][1] = 2;
	node_num[2][2] = 6;
	node_num[2][3] = 5;
	//face 4 high y
	node_num[3][0] = 4;
	node_num[3][1] = 3;
	node_num[3][2] = 7;
	node_num[3][3] = 8;
	//face 5 low z
	node_num[4][0] = 1;
	node_num[4][1] = 2;
	node_num[4][2] = 3;
	node_num[4][3] = 4;
	//face 6 high z
	node_num[5][0] = 5;
	node_num[5][1] = 6;
	node_num[5][2] = 7;
	node_num[5][3] = 8;

	vector< vector<double> > vertices;
	for (int face = 0; face < 6; face++)
	{
	    vertices = mesh->get_vertices(1,face+1);

	    if (vertices.size()    != 3) ITFAILS; // 3 dimensions
	    if (vertices[0].size() != 4) ITFAILS; // 4 nodes per face
	    if (vertices[1].size() != 4) ITFAILS; // 4 nodes per face
	    if (vertices[2].size() != 4) ITFAILS; // 4 nodes per face

	    for (int face_node = 0; face_node < 4; face_node++)
	    {
		for (int dim = 0; dim<3; dim++)
		{
		    if (!soft_equiv(vertices[dim][face_node],
			ref_vert[dim][node_num[face][face_node]-1])) ITFAILS;
		}
	    }
	}
	vertices = mesh->get_vertices(1);

	if (vertices.size()    != 3) ITFAILS; // 3 dimensions
	if (vertices[0].size() != 8) ITFAILS; // 8 nodes per cell
	if (vertices[1].size() != 8) ITFAILS; // 8 nodes per cell
	if (vertices[2].size() != 8) ITFAILS; // 8 nodes per cell

	for (int cell_node = 0; cell_node < 8; cell_node++)
	{
	    for (int dim = 0; dim < 3; dim++)
	    {
		if (!soft_equiv(vertices[dim][cell_node],
				ref_vert[dim][cell_node])) ITFAILS;
	    }
	}
    }

    // checke that the mesh returns the correct face number
    if (mesh->get_bndface("lox",1) != 1) ITFAILS;
    if (mesh->get_bndface("hix",1) != 2) ITFAILS;
    if (mesh->get_bndface("loy",1) != 3) ITFAILS;
    if (mesh->get_bndface("hiy",1) != 4) ITFAILS;
    if (mesh->get_bndface("loz",1) != 5) ITFAILS;
    if (mesh->get_bndface("hiz",1) != 6) ITFAILS;

    // check that, if this were a surface source cell, it would actually be
    // on a vacuum boundary
    {
	vector<int> sur_cell(1,1);
	if (mesh->check_defined_surcells("lox",sur_cell))  ITFAILS;
	if (!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
	if (mesh->check_defined_surcells("loy",sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("hiy",sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("loz",sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("hiz",sur_cell))  ITFAILS;
    }
    {
	vector<int> neighbors;
	neighbors = mesh->get_neighbors(1);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 1) ITFAILS;
	if (neighbors[1]     != 0) ITFAILS;
	if (neighbors[2]     != 1) ITFAILS;
	if (neighbors[3]     != 1) ITFAILS;
	if (neighbors[4]     != 1) ITFAILS;
	if (neighbors[5]     != 1) ITFAILS;
    }

    // make sure the mesh thinks it's not submesh
    if (!mesh->full_Mesh()) ITFAILS;

    // >>>> test the mesh's ability to sample positions. <<<<

    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    if(seed!=1234567)  ITFAILS;
    
    // The explicit-value tests has been hand-checked in a weak sense and,
    // therefore, are more of a regression test.

    // sample positions in cell
    {
	vector<double> position_sampled;
	
	position_sampled = mesh->sample_pos(1,ran_object);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0], 0.7289, 0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1],-0.5642, 0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2], 0.5035, 0.001)) ITFAILS;

	position_sampled = mesh->sample_pos(1,ran_object);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0], 0.6051, 0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1],-0.2852, 0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2],-0.2427, 0.001)) ITFAILS;

	// sample with slope -- explicit-value test;
	vector<double> slope(3);
	slope[0] =-1.0;
	slope[1] = 0.0;
	slope[2] = 0.0;

	position_sampled = mesh->sample_pos(1, ran_object, slope, 1.0);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0],  0.3676,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1], -3.281E-3,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2], -0.2411,0.001)) ITFAILS;

	
	position_sampled = mesh->sample_pos(1, ran_object, slope, 1.0);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0],  0.8785,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1], -0.2575,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2], -0.3002,0.001)) ITFAILS;	
    }

    // sample positions on face
    {
	vector<double> position_sampled;

	position_sampled = mesh->sample_pos_on_face(1,2,ran_object);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0],  1.0   ,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1], -0.5269,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2], -0.0814,0.001)) ITFAILS;

	position_sampled = mesh->sample_pos_on_face(1,2,ran_object);
	if (position_sampled.size() != 3)                    ITFAILS;
	if (!soft_equiv(position_sampled[0],  1.0   ,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[1], -0.8467,0.001)) ITFAILS;
	if (!soft_equiv(position_sampled[2], -0.8006,0.001)) ITFAILS;	
    }

    // now, examine sampled distributions
    // first, sample position using slope
    {
	int num_bins = 100;
	vector<double> xdist(num_bins, 0.0);
	vector<double> ydist(num_bins, 0.0);
	vector<double> zdist(num_bins, 0.0);
	
	double xavg = 0.0;
	double yavg = 0.0;
	double zavg = 0.0;
	double fabs_yavg = 0.0;
	double fabs_zavg = 0.0;

	int xbin, ybin, zbin;
	int num_particles = 100000;

	// -1.0 slope, cell-averaged value of 1 for weighting function
	vector<double> slope(3);
	slope[0] = -1.0;
	slope[1] =  0.0;
	slope[2] =  0.0 ;
	double center_value = 1.0;

	// cell parameters
	double low_y_limit  = -tan_beta;
	double high_y_limit =  tan_beta;

	double xbin_width = 1.0/num_bins;
	double ybin_width = (high_y_limit-low_y_limit)/num_bins;

	// sample position
	vector<double> position_sampled(3);
	for (int num_p = 0; num_p < num_particles; num_p++)
	{
	    position_sampled = mesh->sample_pos(1, ran_object, slope, 
	    					center_value);

	    xbin = static_cast<int>(position_sampled[0]*num_bins);
	    ybin = static_cast<int>((position_sampled[1]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    zbin = static_cast<int>((position_sampled[2]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    xdist[xbin]++;
	    ydist[ybin]++;
	    zdist[zbin]++;

	    xavg += position_sampled[0];
	    yavg += position_sampled[1];
	    zavg += position_sampled[2];

	    fabs_yavg += fabs(position_sampled[1]);
	    fabs_zavg += fabs(position_sampled[2]);
	}

	xavg /= num_particles;
	yavg /= num_particles;
	zavg /= num_particles;

	fabs_yavg /= num_particles;
	fabs_zavg /= num_particles;

	// rough estimate of 1 standard deviation
	double est_std_of_mean = 1.0/sqrt(num_particles);

	// check averages (average value of x is 57/80)
	if (!soft_equiv(xavg, 57./80., 4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(yavg, 0.0    , 4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(zavg, 0.0    , 4.*est_std_of_mean))  ITFAILS;
	
	// check absolute avergages (both fabs_yavg and fabs_zavg
	// = tan_beta*57/160)	
	if (!soft_equiv(fabs_yavg, tan_beta*57/160, 4.*est_std_of_mean))  
	    ITFAILS;
	if (!soft_equiv(fabs_zavg, tan_beta*57/160, 4.*est_std_of_mean))  
	    ITFAILS;

	// now examine distributions themselves
	for (int bin = 0; bin < num_bins; bin++)
	{
	    xdist[bin] /= xbin_width*num_particles;
	    ydist[bin] /= ybin_width*num_particles;
	    zdist[bin] /= ybin_width*num_particles;
	}

	double estimated_std = 
	    sqrt(static_cast<double>(num_bins)/num_particles);

	// check if within 4 sigma
	// x position, ignore first 20 bins (distribution has a tail)
	for (int bin = 20; bin < num_bins; bin++)
	{
	    double xvalue = (bin+0.5)*xbin_width;
	    double fvalue = 21./4.*xvalue*xvalue-3*xvalue*xvalue*xvalue;
	    if (fabs(xdist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	    
	}
	// check y and z distributions
	// ignore first and last 10
	for (int bin = 10; bin < num_bins-10; bin++)
	{
	    double yvalue = low_y_limit+(bin+0.5)*ybin_width;
	    double fvalue = 7./8.*(1.-yvalue*yvalue/(tan_beta*tan_beta));
	    if (yvalue > 0)
	    {
		fvalue -= 1./3.*(1.-yvalue*yvalue*yvalue/
				 (tan_beta*tan_beta*tan_beta));
	    }
	    else
	    {
		fvalue -= 1./3.*(1+yvalue*yvalue*yvalue/
				 (tan_beta*tan_beta*tan_beta));
	    }
	    fvalue *= (3./2.)/tan_beta;

	    if (fabs(ydist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;

	    if (fabs(zdist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	    
	}   
    }
    // next, sample positions using uniform distribution
    {
	int num_bins = 100;
	vector<double> xdist(num_bins, 0.0);
	vector<double> ydist(num_bins, 0.0);
	vector<double> zdist(num_bins, 0.0);

	double xavg = 0.0;
	double yavg = 0.0;
	double zavg = 0.0;
	double fabs_yavg = 0.0;
	double fabs_zavg = 0.0;

	int xbin, ybin, zbin;
	int num_particles = 100000;

	// cell paramters
	double low_y_limit  = -tan_beta;
	double high_y_limit =  tan_beta;

	double xbin_width = 1.0/num_bins;
	double ybin_width = (high_y_limit-low_y_limit)/num_bins;

	// sample position
	vector<double> position_sampled(3);
	for (int num_p = 0; num_p < num_particles; num_p++)
	{
	    position_sampled = mesh->sample_pos(1, ran_object);
 
	    xbin = static_cast<int>(position_sampled[0]*num_bins);
	    ybin = static_cast<int>((position_sampled[1]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    zbin = static_cast<int>((position_sampled[2]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    xdist[xbin]++;
	    ydist[ybin]++;
	    zdist[zbin]++;

	    xavg += position_sampled[0];
	    yavg += position_sampled[1];
	    zavg += position_sampled[2];

	    fabs_yavg += fabs(position_sampled[1]);
	    fabs_zavg += fabs(position_sampled[2]);
	}

	xavg /= num_particles;
	yavg /= num_particles;
	zavg /= num_particles;

	fabs_yavg /= num_particles;
	fabs_zavg /= num_particles;

	// rough estimate of 1 standard deviation
	double est_std_of_mean = 1.0/sqrt(num_particles);

	// check averages (average value of x is 3/4)
	if (!soft_equiv(xavg, 3./4.  , 4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(yavg, 0.0    , 4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(zavg, 0.0    , 4.*est_std_of_mean))  ITFAILS;


	// check absolute avergages (both fabs_yavg and fabs_zavg
	// = 3/8*tan_beta)	
	if (!soft_equiv(fabs_yavg, 3./8.*tan_beta, 
			4.*est_std_of_mean)) ITFAILS;
	if (!soft_equiv(fabs_zavg, 3./8.*tan_beta, 
			4.*est_std_of_mean)) ITFAILS;

	// now examine distribuitions themselves
	for (int bin = 0; bin < num_bins; bin++)
	{
	    xdist[bin] /= xbin_width*num_particles;
	    ydist[bin] /= ybin_width*num_particles;
	    zdist[bin] /= ybin_width*num_particles;
	}

	double estimated_std = 
	    sqrt(static_cast<double>(num_bins)/num_particles);

	// check if within 4 sigma
	// xposition, ignore first 20 bins
	for (int bin = 20; bin < num_bins; bin++)
	{
	    double xvalue = (bin+0.5)*xbin_width;
	    double fvalue =3.*xvalue*xvalue;
	    if (fabs(xdist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	   
	}

	// check y and z distributions
	// ignore first and last ten
	for (int bin = 10; bin < num_bins-10; bin++)
	{
	    double yvalue = low_y_limit+(bin+0.5)*ybin_width;
	    double fvalue = (3./4.)*(1.-yvalue*yvalue/(tan_beta*tan_beta))
		/tan_beta;
	    if (fabs(ydist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	    if (fabs(zdist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	    
	}
    }

    // now test sampling on face

    {
	int num_bins = 100;	
	vector<double> ydist(num_bins, 0.0);
	vector<double> zdist(num_bins, 0.0);
	
	double yavg = 0.0;
	double zavg = 0.0;
	double fabs_yavg = 0.0;
	double fabs_zavg = 0.0;

	int ybin, zbin;
	int num_particles = 100000;
	
	// cell paramters
	double low_y_limit  = -tan_beta;
	double high_y_limit =  tan_beta;

	double ybin_width = (high_y_limit-low_y_limit)/num_bins;

	// sample position
	vector<double> position_sampled(3);
	for (int num_p = 0; num_p < num_particles; num_p++)
	{
	    position_sampled = mesh->sample_pos_on_face(1, 2, ran_object);
 
	    if(!soft_equiv(position_sampled[0], 1.0))  ITFAILS;
	    ybin = static_cast<int>((position_sampled[1]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    zbin = static_cast<int>((position_sampled[2]-low_y_limit)
				    *num_bins/(high_y_limit-low_y_limit));
	    ydist[ybin]++;
	    zdist[zbin]++;

	    yavg += position_sampled[1];
	    zavg += position_sampled[2];

	    fabs_yavg += fabs(position_sampled[1]);
	    fabs_zavg += fabs(position_sampled[2]);
	}

	yavg /= num_particles;
	zavg /= num_particles;

	fabs_yavg /= num_particles;
	fabs_zavg /= num_particles;   

	// rough estimate of 1 standard deviation
	double est_std_of_mean = 1.0/sqrt(num_particles);

	// check averages (average value of fabs_y and fabs_z is 0.5*tan_beta)
	if (!soft_equiv(yavg,      0.0,          4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(zavg,      0.0,          4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(fabs_yavg, 0.5*tan_beta, 4.*est_std_of_mean))  ITFAILS;
	if (!soft_equiv(fabs_zavg, 0.5*tan_beta, 4.*est_std_of_mean))  ITFAILS;

	// now examine distribuitions themselves
	for (int bin = 0; bin < num_bins; bin++)
	{
	    ydist[bin] /= ybin_width*num_particles;
	    zdist[bin] /= ybin_width*num_particles;
	}
	
	double estimated_std = 
	    sqrt(static_cast<double>(num_bins)/num_particles);

	// check if within 4 sigma
	// check y and z distributions
	// ignore first and last ten
	for (int bin = 0; bin < num_bins-0; bin++)
	{
	    double fvalue = 0.5/tan_beta;
	    if (fabs(ydist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
 	    if (fabs(zdist[bin]-fvalue) > 4.*estimated_std*fvalue) ITFAILS;
	}
    }

    // check the == and != operations.
    // first, build another mesh object equivalent to old mesh object.
    // these two meshes, althoug identical, should not occupy the same memory
    {
	SP<Sphyramid_Mesh> other_mesh(new Sphyramid_Mesh(coord, layout,
							 cell_x_extents,
							 beta_radians));
	if (mesh == other_mesh) ITFAILS;

	Sphyramid_Mesh &mesh_obj       = *mesh;
	Sphyramid_Mesh &other_mesh_obj = *other_mesh;

	if (mesh_obj   != other_mesh_obj)  ITFAILS;
	if (!(mesh_obj == other_mesh_obj)) ITFAILS;
    }

    // test the mesh's distance-to-boundary function
    {
	vector<double> r(3,0.0);
	vector<double> omega(3,0.0);
	int intersecting_face;
	double db;

	omega[0] = 1.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 1.0))   ITFAILS;
	if (intersecting_face != 2) ITFAILS;

	r[0]     =  0.5;
	r[1]     = -0.5*tan_beta;
	r[2]     =  0.2;
	omega[0] =  0.0;
	omega[1] =  1.0;
	omega[2] =  0.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, tan_beta)) ITFAILS;
	if (intersecting_face != 4)    ITFAILS;

	r[0]     =  0.5;
	r[1]     =  0.0;
	r[2]     =  0.2;
	omega[0] =  0.0;
	omega[1] =  0.0;
	omega[2] =  1.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 0.5*tan_beta-0.2)) ITFAILS;
	if (intersecting_face != 6)            ITFAILS;

	r[0]     =  0.5;
	r[1]     =  0.0;
	r[2]     =  0.0;
	omega[0] = -0.70710678119;
	omega[1] =  0.70710678119;
	omega[2] =  0.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 0.332227824872, 1.E-10)) ITFAILS;
	if (intersecting_face != 4)                  ITFAILS;

	r[0]     =  0.5;
	r[1]     =  0.0;
	r[2]     =  0.0;
	omega[0] = -0.70710678119;
	omega[1] = -0.70710678119;
	omega[2] =  0.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 0.332227824872, 1.E-10)) ITFAILS;
	if (intersecting_face != 3)                  ITFAILS;
	
	r[0]     =  0.5678;
	r[1]     =  0.0;
	r[2]     =  0.0;
	omega[0] = -1.0;
	omega[1] =  0.0;
	omega[2] =  0.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 0.5678, 1.E-10)) ITFAILS;
	if (!((intersecting_face == 1) || (intersecting_face == 3) ||
	      (intersecting_face == 4) || (intersecting_face == 5) ||
	      (intersecting_face == 6))) ITFAILS;

	r[0]     =  1.0;
	r[1]     = -0.5;
	r[2]     =  0.1;
	omega[0] = -1.0;
	omega[1] =  0.0;
	omega[2] =  0.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, 1.0-0.5/tan_beta, 1.E-10)) ITFAILS;
	if (intersecting_face != 3)                    ITFAILS;
	
	r[0]     =  0.5;
	r[1]     =  0.0;
	r[2]     =  0.5*tan_beta;
	omega[0] =  0.0;
	omega[1] =  0.0;
	omega[2] = -1.0;
	db = mesh->get_db(r, omega, 1, intersecting_face);
	if (!soft_equiv(db, tan_beta, 1.E-10)) ITFAILS;
	if (intersecting_face != 5)            ITFAILS;

    }

    // test get_random_walk_sphere_radius fucntion
    {
	// starting point
	vector<double> r(3);
	r[0] = 0.5;
	r[1] = 0.0;
	r[2] = 0.0;
	if (!mesh->in_cell(1, r)) ITFAILS;

	// get random walk sphere radius for this position
	double rw_radius = mesh->get_random_walk_sphere_radius(r, 1);

	if (!soft_equiv(rw_radius, 0.5-1.0E-6)) ITFAILS;

	// move point and try again
	r[0] -= 0.25;

	rw_radius = mesh->get_random_walk_sphere_radius(r, 1);

	if (!soft_equiv(rw_radius, 0.25-1.0E-6)) ITFAILS;
    }

    // test sample position on sphere
    {
	// make a random number generator
	Rnd_Control control(seed);
	Sprng ran     = control.get_rn(10);
	Sprng ref_ran = control.get_rn(10);

	// random walk sphere radius
	double radius;

	// starting position
        vector<double> r(3);
	r[0] = 0.5;
	r[1] = 0.0;
	r[2] = 0.0;

	// new and reference position and direction
	vector<double> r_prime(3);
	vector<double> omega_prime(3);
	vector<double> ref_r(3);
	vector<double> ref_omega(3);

	pair<vector<double>, vector<double> > r_and_omega;

	// simulate a bunch of random walk steps
	for (int j = 0; j < 20; j++)
	{

	    // save a reference position
	    ref_r = r;

	    // calculate the RW radius for this starting point
	    radius = mesh->get_random_walk_sphere_radius(r, 1);

	    // sample new position and direction
	    r_and_omega = mesh->sample_random_walk_sphere(1, r, radius, ran);
	    r_prime     = r_and_omega.first;
	    omega_prime = r_and_omega.second;
	    
	    // check that new position is in cell
	    if (!mesh->in_cell(1, r_prime)) ITFAILS;

	    // calculate reference isotropic direction (equivalent to the
	    // first random walk direction each time)
	    ref_omega = mesh->get_Coord().sample_isotropic_dir(ref_ran);

	    // if the angle are the same there were no y or z face directions
	    if(soft_equiv(omega_prime.begin(), omega_prime.end(),
			  ref_omega.begin(), ref_omega.end()))
	    {
		// check that sphere is inside cell
		int face;
		double db;
		db = mesh->get_db(r, ref_omega, 1, face);
		if (radius > db) ITFAILS;

		// update position
		ref_r[0]+=ref_omega[0]*radius;
		ref_r[1]+=ref_omega[1]*radius;
		ref_r[2]+=ref_omega[2]*radius;

		if (!soft_equiv(r_prime.begin(), r_prime.end(),
				ref_r.begin(), ref_r.end())) ITFAILS;
	    }

	    // else there were y or z reflections
	    else
	    {
		// check that y or z face is hit first

		int face;
		double db;
		db = mesh->get_db(r, ref_omega, 1, face);
		if (face   < 3)   ITFAILS;
		if (face   > 6)   ITFAILS;
		if (radius < db ) ITFAILS;
	    }

	    // check that the net distance traveled <= radius
	    double dist = sqrt( (r_prime[0]-r[0])*(r_prime[0]-r[0])
				+ (r_prime[1]-r[1])*(r_prime[1]-r[1])
				+ (r_prime[2]-r[2])*(r_prime[2]-r[2]));

	    if (!soft_equiv(dist, radius))
	    { 
		if (dist > radius) ITFAILS;
	    }

	    // reassign postion;
	    
	    r = r_prime;
	}
    }

    return;
}
//---------------------------------------------------------------------------//
// test the Sphyramid_Mesh via the Sphyramid_Builder
void build_a_Sphyramid()
{
    using rtt_dsxx::SP;
    using rtt_mc_test::Parser;
    using rtt_mc::Sphyramid_Builder;
    using rtt_mc::Sphyramid_Mesh;
    using std::vector;
    using std::string;
    using rtt_mc::global::pi;
    using std::pow;
    using rtt_dsxx::soft_equiv;
    using std::cos;
    using std::tan;
    using std::sin;
    using std::atan;

    SP<Parser> parser(new Parser("Sphyramid_Input"));
    Sphyramid_Builder builder(parser);

    // check some of the Sphyramid_Builder properties.

    // test cell region data
    {
	vector<int> regions(5,1);
	regions[2] = 2;
	regions[3] = 2;
	regions[4] = 2;

	if (builder.get_num_regions() != 2)       ITFAILS;
	if (builder.get_regions()     != regions) ITFAILS;
    }

    // test zone mapping
    {
	vector< vector<int> > zone;
	zone.resize(2);
	zone[0].resize(2);
	zone[1].resize(3);

	zone[0][0] = 1;
	zone[0][1] = 2;
	zone[1][0] = 3;
	zone[1][1] = 4;
	zone[1][2] = 5;

	if (builder.get_num_zones()      != 2)       ITFAILS;
	if (builder.get_cells_in_zone(1) != zone[0]) ITFAILS;
	if (builder.get_cells_in_zone(2) != zone[1]) ITFAILS;

    }
    
    //build a Sphyramid_mesh object
    SP<Sphyramid_Mesh> mesh = builder.build_Mesh();

    //check defined surface source cells
    {
	vector< vector<int> > ss(1);
	ss[0].resize(1);
	ss[0][0] = 5;

	if (builder.get_defined_surcells() != ss) ITFAILS;

	vector<string> ssp(1);
	ssp[0]="hir";
	if (builder.get_ss_pos() != ssp) ITFAILS;
    }

    // check zone mapper
    {
	vector<int> zone_field(2);
	zone_field[0] = 1000;
	zone_field[1] = 1001;
	vector<int> cell_field = builder.zone_cell_mapper(zone_field);

	if (cell_field.size() != mesh->num_cells()) ITFAILS;
	if (cell_field[0]     != 1000)              ITFAILS;
	if (cell_field[1]     != 1000)              ITFAILS;
	if (cell_field[2]     != 1001)              ITFAILS;
	if (cell_field[3]     != 1001)              ITFAILS;
	if (cell_field[4]     != 1001)              ITFAILS;
    }

    // check number of cells
    if (builder.num_cells() != 5) ITFAILS;

    // check get_Mesh
    {
	SP<Sphyramid_Mesh> same_mesh = builder.get_Mesh();
	if (same_mesh  !=  mesh) ITFAILS;
	if (*same_mesh != *mesh) ITFAILS;
    }

    // <<<< now do some checks on the mesh itself >>>>
    // (this is a different mesh than the simple_one_cell_Sphyramid_Mesh
    
    // check that the mesh returns proper coordinate system information
    if (mesh->get_Coord().get_dim()      != 3)             ITFAILS;
    if (mesh->get_SPCoord()->get_dim()   != 3)             ITFAILS;
    if (mesh->get_Coord().get_Coord()    != string("xyz")) ITFAILS;
    if (mesh->get_SPCoord()->get_Coord() != string("xyz")) ITFAILS;

    // check that the mesh returns the proper number of cells
    if (mesh->num_cells() != 5) ITFAILS;

    // calculate mesh properties that will be needed throughout the remainder
    // of the unit test

    // spherical cone angle (alpha) is 45 degrees
    double r_to_x   = pow((2.*(1.-cos(pi/4.)))/(tan(pi/4)*tan(pi/4)), 1./3.);
    double tan_beta = sqrt(pi)/2.;
    double cos_beta = cos(atan(tan_beta));	
    double sin_beta = sin(atan(tan_beta));

    // the radial boundaries are 0,0.5,1,9/7,13/7,3
    vector<double> xcoords(6);
    xcoords[0] = 0.0*r_to_x;
    xcoords[1] = 0.5*r_to_x;
    xcoords[2] = 1.0*r_to_x;
    xcoords[3] = 9./7.*r_to_x;
    xcoords[4] = 13./7.*r_to_x;
    xcoords[5] = 3.*r_to_x;

    // check that the x extents are correct for the 1st and last cells
    if (!soft_equiv(mesh->get_low_x(1),  xcoords[0])) ITFAILS;
    if (!soft_equiv(mesh->get_high_x(1), xcoords[1])) ITFAILS; 

    if (!soft_equiv(mesh->get_low_x(5),  xcoords[4])) ITFAILS;
    if (!soft_equiv(mesh->get_high_x(5), xcoords[5])) ITFAILS;

    // check that cell midpoints are correct for 1st and last cell
    {
	double x_midpoint = 0.5*(xcoords[0]+xcoords[1]);
	if (!soft_equiv(mesh->get_x_midpoint(1), x_midpoint)) ITFAILS;
	if (!soft_equiv(mesh->get_y_midpoint(1), 0.0))        ITFAILS;
	if (!soft_equiv(mesh->get_z_midpoint(1), 0.0))        ITFAILS;
    }
    {
	double x_midpoint = 0.5*(xcoords[4]+xcoords[5]);
	if (!soft_equiv(mesh->get_x_midpoint(5), x_midpoint)) ITFAILS;
	if (!soft_equiv(mesh->get_y_midpoint(5), 0.0))        ITFAILS;
	if (!soft_equiv(mesh->get_z_midpoint(5), 0.0))        ITFAILS;
    }
    
    // check that the cell dimensions are correct for 1st and last cells;
    {
	double x_dim = xcoords[1]-xcoords[0];
	double y_dim = 0.5*(xcoords[0]+xcoords[1])*2.*tan_beta;
	double z_dim = y_dim;

	if (!soft_equiv(mesh->dim(1,1), x_dim)) ITFAILS;
	if (!soft_equiv(mesh->dim(2,1), y_dim)) ITFAILS;
	if (!soft_equiv(mesh->dim(3,1), z_dim)) ITFAILS;
    }
    {
	double x_dim = xcoords[5]-xcoords[4];
	double y_dim = 0.5*(xcoords[5]+xcoords[4])*2.*tan_beta;
	double z_dim = y_dim;

	if (!soft_equiv(mesh->dim(1,5), x_dim)) ITFAILS;
	if (!soft_equiv(mesh->dim(2,5), y_dim)) ITFAILS;
	if (!soft_equiv(mesh->dim(3,5), z_dim)) ITFAILS;
    }

    // check that graphics cell type is correct
    {
	vector<int> cell_type = mesh->get_cell_types();
	if (cell_type.size() != 5)                              ITFAILS;
	for (int i = 0; i < 5; i++)
	{
	    if (cell_type[0] != rtt_viz::eight_node_hexahedron) ITFAILS;
	}
    }

    // check that mesh returns the proper point coordinates (1st & last cell)
    {
	vector< vector<double> > point_coords = mesh->get_point_coord();
	if (point_coords.size() != 40)                             ITFAILS;
	for (int p = 0; p < 40; p++)
	{
	    if (point_coords[p].size() != 3)                       ITFAILS;
	}

	// first cell, low z face, counter clockwise from the lower left, 
	// facing +z
	if (!soft_equiv(point_coords[0][0], xcoords[0]))           ITFAILS;
	if (!soft_equiv(point_coords[0][1],-xcoords[0]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[0][2],-xcoords[0]*tan_beta))  ITFAILS;

	if (!soft_equiv(point_coords[1][0], xcoords[1]))           ITFAILS;
	if (!soft_equiv(point_coords[1][1],-xcoords[1]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[1][2],-xcoords[1]*tan_beta))  ITFAILS;

	if (!soft_equiv(point_coords[2][0], xcoords[1]))           ITFAILS;
	if (!soft_equiv(point_coords[2][1], xcoords[1]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[2][2],-xcoords[1]*tan_beta))  ITFAILS; 

	if (!soft_equiv(point_coords[3][0], xcoords[0]))           ITFAILS;
	if (!soft_equiv(point_coords[3][1], xcoords[0]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[3][2],-xcoords[0]*tan_beta))  ITFAILS;

	// first cell, high z face, counter clocwise from the lower left, 
	// facing +z
	if (!soft_equiv(point_coords[4][0], xcoords[0]))           ITFAILS;
	if (!soft_equiv(point_coords[4][1],-xcoords[0]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[4][2], xcoords[0]*tan_beta))  ITFAILS;

	if (!soft_equiv(point_coords[5][0], xcoords[1]))           ITFAILS;
	if (!soft_equiv(point_coords[5][1],-xcoords[1]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[5][2], xcoords[1]*tan_beta))  ITFAILS;

	if (!soft_equiv(point_coords[6][0], xcoords[1]))           ITFAILS;
	if (!soft_equiv(point_coords[6][1], xcoords[1]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[6][2], xcoords[1]*tan_beta))  ITFAILS;

	if (!soft_equiv(point_coords[7][0], xcoords[0]))           ITFAILS;
	if (!soft_equiv(point_coords[7][1], xcoords[0]*tan_beta))  ITFAILS;
	if (!soft_equiv(point_coords[7][2], xcoords[0]*tan_beta))  ITFAILS;

	// last cell, low z face, counter clockwise from the lower left, 
	// facing +z	
	if (!soft_equiv(point_coords[32][0], xcoords[4]))          ITFAILS;
	if (!soft_equiv(point_coords[32][1],-xcoords[4]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[32][2],-xcoords[4]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[33][0], xcoords[5]))          ITFAILS;	
	if (!soft_equiv(point_coords[33][1],-xcoords[5]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[33][2],-xcoords[5]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[34][0], xcoords[5]))          ITFAILS;
	if (!soft_equiv(point_coords[34][1], xcoords[5]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[34][2],-xcoords[5]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[35][0], xcoords[4]))          ITFAILS;
	if (!soft_equiv(point_coords[35][1], xcoords[4]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[35][2],-xcoords[4]*tan_beta)) ITFAILS;
	
	// last cell, high z face, counter clockwise from the lower left, 
	// facing +z
	if (!soft_equiv(point_coords[36][0], xcoords[4]))          ITFAILS;
	if (!soft_equiv(point_coords[36][1],-xcoords[4]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[36][2], xcoords[4]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[37][0], xcoords[5]))          ITFAILS;
	if (!soft_equiv(point_coords[37][1],-xcoords[5]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[37][2], xcoords[5]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[38][0], xcoords[5]))          ITFAILS;
	if (!soft_equiv(point_coords[38][1], xcoords[5]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[38][2], xcoords[5]*tan_beta)) ITFAILS;

	if (!soft_equiv(point_coords[39][0], xcoords[4]))          ITFAILS;
	if (!soft_equiv(point_coords[39][1], xcoords[4]*tan_beta)) ITFAILS;
	if (!soft_equiv(point_coords[39][2], xcoords[4]*tan_beta)) ITFAILS;
    }

    // check that mesh returns the correct "next cell"
    // low and high y and z should be reflecting
    for (int c = 1; c <= 5; c++)
    {
	if (mesh->next_cell(c,3) != c) ITFAILS;
	if (mesh->next_cell(c,4) != c) ITFAILS;
	if (mesh->next_cell(c,5) != c) ITFAILS;
	if (mesh->next_cell(c,6) != c) ITFAILS;
    }
    if (mesh->next_cell(1,1) != 1)     ITFAILS;
    if (mesh->next_cell(1,2) != 2)     ITFAILS;
    if (mesh->next_cell(2,1) != 1)     ITFAILS;
    if (mesh->next_cell(2,2) != 3)     ITFAILS;
    if (mesh->next_cell(3,1) != 2)     ITFAILS;
    if (mesh->next_cell(3,2) != 4)     ITFAILS;
    if (mesh->next_cell(4,1) != 3)     ITFAILS;
    if (mesh->next_cell(4,2) != 5)     ITFAILS;
    if (mesh->next_cell(5,1) != 4)     ITFAILS;
    if (mesh->next_cell(5,2) != 0)     ITFAILS;
    
    // check that mesh can locate a cell given a position
    {
	vector<double> position(3,0.0);
	position[0] = 0.25*r_to_x;
	if (mesh->get_cell(position) != 1) ITFAILS;
   
	position[0] = 100.0;
	if (mesh->get_cell(position) !=-1) ITFAILS;

	position[0] = 10./7.*r_to_x;
	position[1] = tan_beta;
	position[2] =-tan_beta;
	if (mesh->get_cell(position) != 4) ITFAILS;
    }

    // check that the face normals are correct
    {
	vector<double> normal(3,0.0);

	normal = mesh->get_normal(1,1);
	if (!soft_equiv(normal[0],-1.0))       ITFAILS;
	if (!soft_equiv(normal[1], 0.0))       ITFAILS;
	if (!soft_equiv(normal[2], 0.0))       ITFAILS;
    
	normal=mesh->get_normal(1,2);
	if (!soft_equiv(normal[0], 1.0))       ITFAILS;
	if (!soft_equiv(normal[1], 0.0))       ITFAILS;
	if (!soft_equiv(normal[2], 0.0))       ITFAILS;
    
	normal=mesh->get_normal(1,3);
	if (!soft_equiv(normal[0],-sin_beta))  ITFAILS;
	if (!soft_equiv(normal[1],-cos_beta))  ITFAILS;
	if (!soft_equiv(normal[2], 0.0))       ITFAILS;
	
	normal=mesh->get_normal(1,4);
	if (!soft_equiv(normal[0],-sin_beta))  ITFAILS;
        if (!soft_equiv(normal[1], cos_beta))  ITFAILS;
        if (!soft_equiv(normal[2], 0.0))       ITFAILS;
	
	normal=mesh->get_normal(1,5);
	if (!soft_equiv(normal[0],-sin_beta))  ITFAILS;
	if (!soft_equiv(normal[1], 0.0))       ITFAILS;
	if (!soft_equiv(normal[2],-cos_beta))  ITFAILS;
	
	normal=mesh->get_normal(1,6);
	if (!soft_equiv(normal[0],-sin_beta))  ITFAILS;
	if (!soft_equiv(normal[1], 0.0))       ITFAILS;
	if (!soft_equiv(normal[2], cos_beta))  ITFAILS;

	// now see that all other cells have the same normals
	for (int c = 1; c <= 5; c++)
	{
	    for (int f = 1; f <= 6; f++)
	    {
		normal = mesh->get_normal(1,f);
		vector<double> check = mesh->get_normal(c,f);
		for (int d = 0; d < 3;d++)
		{
		    if (!soft_equiv(normal[d], check[d]))  ITFAILS;
		}
	    }
	}

	vector<double> normal_in(3,0.0);
	for (int cell = 1; cell <= 5; cell++)
	{
	    for (int face = 1; face <= 6; face++)
	    {
		normal    = mesh->get_normal(1,face);
		normal_in = mesh->get_normal_in(1,face);

		for (int dim = 0; dim < 3; dim++)
		{
		    if (!soft_equiv(normal[dim],-normal_in[dim])) ITFAILS;
		}
	    }
	}
    }

    // does the mesh return the correct cell volume and total volume?
    {    
	double cell_volume;
	double total_volume;
	const double vol_factor = 4./3.*tan_beta*tan_beta;
	
	cell_volume  = vol_factor*(pow(xcoords[1],3)-pow(xcoords[0],3));
	if (!soft_equiv(mesh->volume(1), cell_volume))           ITFAILS;

	cell_volume  = vol_factor*(pow(xcoords[2],3)-pow(xcoords[1],3));
	if (!soft_equiv(mesh->volume(2), cell_volume))           ITFAILS;

	cell_volume  = vol_factor*(pow(xcoords[3],3)-pow(xcoords[2],3));
	if (!soft_equiv(mesh->volume(3), cell_volume))           ITFAILS;

	cell_volume  = vol_factor*(pow(xcoords[4],3)-pow(xcoords[3],3));
	if (!soft_equiv(mesh->volume(4), cell_volume))           ITFAILS;

	cell_volume  = vol_factor*(pow(xcoords[5],3)-pow(xcoords[4],3));
	if (!soft_equiv(mesh->volume(5), cell_volume))           ITFAILS;

	total_volume = vol_factor*(pow(xcoords[5],3)-pow(xcoords[0],3));
	if (!soft_equiv(mesh->get_total_volume(), total_volume)) ITFAILS;
    }

    // does the mesh return the correct face area?
    {
	double area;

	// cell 1
	area = pow(2.*xcoords[0]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(1,1), area)) ITFAILS;
    
	area = pow(2.*xcoords[1]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(1,2), area)) ITFAILS;

	area  = xcoords[1]*xcoords[1]-xcoords[0]*xcoords[0];
	area *= tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(1,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(1,6), area)) ITFAILS;

	// cell 2
	area = pow(2.*xcoords[1]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(2,1), area)) ITFAILS;
    
	area = pow(2.*xcoords[2]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(2,2), area)) ITFAILS;

	area  = xcoords[2]*xcoords[2]-xcoords[1]*xcoords[1];
	area *= tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(2,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(2,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(2,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(2,6), area)) ITFAILS;

	// cell 3
	area = pow(2.*xcoords[2]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(3,1), area)) ITFAILS;
    
	area = pow(2.*xcoords[3]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(3,2), area)) ITFAILS;

	area  = xcoords[3]*xcoords[3]-xcoords[2]*xcoords[2];
	area *= tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(3,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(3,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(3,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(3,6), area)) ITFAILS;
	
	// cell 4
	area = pow(2.*xcoords[3]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(4,1), area)) ITFAILS;
    
	area = pow(2.*xcoords[4]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(4,2), area)) ITFAILS;

	area  = xcoords[4]*xcoords[4]-xcoords[3]*xcoords[3];
	area *= tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(4,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(4,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(4,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(4,6), area)) ITFAILS;

	// cell 5
	area = pow(2.*xcoords[4]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(5,1), area)) ITFAILS;
    
	area = pow(2.*xcoords[5]*tan_beta,2);
	if (!soft_equiv(mesh->face_area(5,2), area)) ITFAILS;

	area  = xcoords[5]*xcoords[5]-xcoords[4]*xcoords[4];
	area *= tan_beta/cos_beta;
	if (!soft_equiv(mesh->face_area(5,3), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(5,4), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(5,5), area)) ITFAILS;
	if (!soft_equiv(mesh->face_area(5,6), area)) ITFAILS;
    }
    
    // check that the mesh returns correctly sized vertices
    {
	vector< vector<double> > vertices;
	for (int cell = 0; cell < 5; cell++)
	{
	    for (int face = 0 ; face < 6; face++)
	    {
		vertices = mesh->get_vertices(cell+1,face+1);

		if (vertices.size()    != 3 ) ITFAILS; // 3 dimensions
		if (vertices[0].size() != 4)  ITFAILS; // 4 nodes per face
		if (vertices[1].size() != 4)  ITFAILS; // 4 nodes per face
		if (vertices[2].size() != 4)  ITFAILS; // 4 nodes per face
	    }

	    vertices = mesh->get_vertices(cell+1);
	    if (vertices.size()    != 3) ITFAILS;  // 3 dimensions
	    if (vertices[0].size() != 8) ITFAILS;  // 8 nodes per cell
	    if (vertices[1].size() != 8) ITFAILS;  // 8 nodes per cell
	    if (vertices[2].size() != 8) ITFAILS;  // 8 nodes per cell

	}

	// spot check a few vertices
	vertices = mesh->get_vertices(5,5);
	if (!soft_equiv(vertices[0][0], xcoords[4]))          ITFAILS;
	if (!soft_equiv(vertices[1][0],-xcoords[4]*tan_beta)) ITFAILS;
	if (!soft_equiv(vertices[2][0],-xcoords[4]*tan_beta)) ITFAILS;
	
	vertices = mesh->get_vertices(2,2);
	if (!soft_equiv(vertices[0][2], xcoords[2]))          ITFAILS;
	if (!soft_equiv(vertices[1][2], xcoords[2]*tan_beta)) ITFAILS;
	if (!soft_equiv(vertices[1][2], xcoords[2]*tan_beta)) ITFAILS;
    }

    // check that the mesh returns the correct face number
    if (mesh->get_bndface("lox", 1) != 1) ITFAILS;
    if (mesh->get_bndface("hix", 1) != 2) ITFAILS;
    if (mesh->get_bndface("lor", 1) != 1) ITFAILS;
    if (mesh->get_bndface("hir", 1) != 2) ITFAILS;
    if (mesh->get_bndface("loy", 1) != 3) ITFAILS;
    if (mesh->get_bndface("hiy", 1) != 4) ITFAILS;
    if (mesh->get_bndface("loz", 1) != 5) ITFAILS;
    if (mesh->get_bndface("hiz", 1) != 6) ITFAILS;

    // check that surface source cells would actuall be on a vacuum boundary
    {
	vector<int> sur_cell;
	sur_cell.resize(1);
	sur_cell[0] = 5;
	if (!mesh->check_defined_surcells("hix", sur_cell)) ITFAILS;
	if (!mesh->check_defined_surcells("hir", sur_cell)) ITFAILS;
    
	sur_cell[0] = 1;
	if (mesh->check_defined_surcells("lox", sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("lor", sur_cell))  ITFAILS;
    
	sur_cell.resize(5);
	for(int c = 0; c < 5; c++)
	{
	    sur_cell[c] = c+1;
	}
	if (mesh->check_defined_surcells("loy", sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("hiy", sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("loz", sur_cell))  ITFAILS;
	if (mesh->check_defined_surcells("hiz", sur_cell))  ITFAILS;
    }

    // check that the mesh returns the proper cell list
    {
	vector<int> surface_cell_list;
	surface_cell_list = mesh->get_surcells("hir");
	if (surface_cell_list.size() != 1) ITFAILS;
	if (surface_cell_list[0]     != 5) ITFAILS;
    }

    // check that the neighbors are correct
    {
	vector<int> neighbors;
	neighbors = mesh->get_neighbors(1);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 1) ITFAILS;
	if (neighbors[1]     != 2) ITFAILS;
	if (neighbors[2]     != 1) ITFAILS;
	if (neighbors[3]     != 1) ITFAILS;
	if (neighbors[4]     != 1) ITFAILS;
	if (neighbors[5]     != 1) ITFAILS;
    
	neighbors = mesh->get_neighbors(2);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 1) ITFAILS;
	if (neighbors[1]     != 3) ITFAILS;
	if (neighbors[2]     != 2) ITFAILS;
	if (neighbors[3]     != 2) ITFAILS;
	if (neighbors[4]     != 2) ITFAILS;
	if (neighbors[5]     != 2) ITFAILS;
    
	neighbors = mesh->get_neighbors(3);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 2) ITFAILS;
	if (neighbors[1]     != 4) ITFAILS;
	if (neighbors[2]     != 3) ITFAILS;
	if (neighbors[3]     != 3) ITFAILS;
	if (neighbors[4]     != 3) ITFAILS;
	if (neighbors[5]     != 3) ITFAILS;
    
	neighbors = mesh->get_neighbors(4);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 3) ITFAILS;
	if (neighbors[1]     != 5) ITFAILS;
	if (neighbors[2]     != 4) ITFAILS;
	if (neighbors[3]     != 4) ITFAILS;
	if (neighbors[4]     != 4) ITFAILS;
	if (neighbors[5]     != 4) ITFAILS;
    
	neighbors = mesh->get_neighbors(5);
	if (neighbors.size() != 6) ITFAILS;
	if (neighbors[0]     != 4) ITFAILS;
	if (neighbors[1]     != 0) ITFAILS;
	if (neighbors[2]     != 5) ITFAILS;
	if (neighbors[3]     != 5) ITFAILS;
	if (neighbors[4]     != 5) ITFAILS;
	if (neighbors[5]     != 5) ITFAILS;
    }

    // make sure the mesh thinks it's not a submesh
    if (!mesh->full_Mesh()) ITFAILS;

    // check in cell function
    // cell 3
    {
	vector<double> r(3);
	r[0] = xcoords[3]-1.E-13;
	r[1] = tan_beta;
	r[2] = tan_beta;

	if (!mesh->in_cell(3,r))  ITFAILS;
    }

    return;
}
//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    using C4::Init;
    using C4::Finalize;
    using C4::node;
    using rtt_mc::release;
    using rtt_dsxx::assertion;
    using rtt_mc_test::passed;
    using std::string;
    using std::cout;
    using std::endl;

    Init(argc, argv);

    // this is a serial test
    if (node())
    {
	Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (node() == 0)
		cout << argv[0] << ": version " << release() 
		     << endl;
	    Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	
	// simple one-celled test problem
	simple_one_cell_Sphyramid();

	// test the Sphyramid mesh Builder
	build_a_Sphyramid();
    }
    catch (assertion &ass)
    {
	cout << "While testing tstSphyramid_Mesh, " << ass.what()
	     << endl;
	Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (passed) 
	{
	    cout << "**** tstSphyramid_Mesh Test: PASSED on " 
		 << node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    cout << "Done testing tstSphyramid_Mesh on " << node() << endl;
    
    Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstSphyramid_Mesh.cc
//---------------------------------------------------------------------------//
