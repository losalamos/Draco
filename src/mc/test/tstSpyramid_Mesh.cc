//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSPyramid_Mesh.cc
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
#include "../Pyramid_Mesh.hh"
#include "../XYZCoord_sys.hh"
#include "../AMR_Layout.hh"
#include "../Spyramid_Builder.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "viz/Ensight_Translator.hh"
#include "rng/Rnd_Control.hh"
#include "c4/global.hh"
//#include "c4/SpinLock.hh"
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
void simple_one_cell_Pyramid()
{
    using rtt_dsxx::SP;
    using rtt_mc::XYZCoord_sys;
    using rtt_mc::AMR_Layout;
    using rtt_mc::Pyramid_Mesh;
    using std::string;
    using std::vector;
    using rtt_dsxx::soft_equiv;
    using std::tan;
    using std::sin;
    using std::cos;
    using rtt_rng::Rnd_Control;
    using rtt_rng::Sprng;


    // >>> one-cell problem
    int ncells =1 ;
    
    //>>> build an XYZ coordinate system class <<<
    SP<XYZCoord_sys> coord(new XYZCoord_sys());

    // >>> build a layout of size ncells <<<
    AMR_Layout layout;
    layout.set_size(ncells);

    // set number of faces per cell (6)
    for (int cell=1; cell <= ncells; cell++)
	layout.set_size(cell,6);

    // set vacuum boundary condition on external boundaries;
    // reflecting on inherently reflecting boundaries
    for(int cell=1; cell <= ncells; cell++)
    {
	layout(cell,1,1)=cell;
	layout(cell,2,1)=0;
	layout(cell,3,1)=cell;
	layout(cell,4,1)=cell;
	layout(cell,5,1)=cell;
	layout(cell,6,1)=cell;
    }

    // >>> set x cell extents <<<
    vector<vector<double> > cell_x_extents;
    cell_x_extents.resize(ncells);
    for(int cell=1;cell <= ncells; cell++)
	cell_x_extents[cell-1].resize(2);

    for(int cell=1; cell <=ncells; cell++)
    {
	cell_x_extents[cell-1][0]=0.0; // low x
	cell_x_extents[cell-1][1]=1.0; //high x
    }

    // set pyramid_angle
    double beta_radians=0.725153345774; // alpha_degrees=45
    
    // calculate tan(beta), you'll need it
    double tan_beta = tan(beta_radians);

    // build a mesh object
    SP<Pyramid_Mesh> mesh(new Pyramid_Mesh(coord,layout,cell_x_extents, 
					   beta_radians));
    
    // check that the mesh returns proper coordinate system information
    if(mesh->get_Coord().get_dim() !=3) ITFAILS;
    if(mesh->get_SPCoord()->get_dim() !=3) ITFAILS;
    if(mesh->get_Coord().get_Coord()!=string("xyz")) ITFAILS;
    if(mesh->get_SPCoord()->get_Coord() !=string("xyz")) ITFAILS;

    // check that mesh returns the proper number of cells
    if(mesh->num_cells() !=ncells)  ITFAILS;

    // check that x extents are correct
    if(!soft_equiv(mesh->get_low_x(1), 0.0))  ITFAILS;
    if(!soft_equiv(mesh->get_high_x(1),1.0))  ITFAILS;

    // check that cell midpoints are correct
    if (!soft_equiv(mesh->get_x_midpoint(1),0.5)) ITFAILS;
    if (!soft_equiv(mesh->get_y_midpoint(1),0.0)) ITFAILS;
    if (!soft_equiv(mesh->get_z_midpoint(1),0.0)) ITFAILS;

    // check that the cell dimensions are correct
    if(!soft_equiv(mesh->dim(1,1), 1.0)) ITFAILS;
    if(!soft_equiv(mesh->dim(2,1), tan_beta )) ITFAILS;
    if(!soft_equiv(mesh->dim(3,1), tan_beta )) ITFAILS;

    // check that graphics cell type is correct
    vector<int> cell_type = mesh->get_cell_types();
    if (cell_type.size()!=1)                            ITFAILS;
    if (cell_type[0] !=rtt_viz::eight_node_hexahedron ) ITFAILS;

    // check that the mesh returns the proper point coordinates
    vector<vector<double> > point_coords = mesh->get_point_coord();
    if(point_coords.size() !=8 )  ITFAILS;
    if(point_coords[0].size()!=3) ITFAILS;
    if(point_coords[1].size()!=3) ITFAILS;
    if(point_coords[2].size()!=3) ITFAILS;
    if(point_coords[3].size()!=3) ITFAILS;
    if(point_coords[4].size()!=3) ITFAILS;
    if(point_coords[5].size()!=3) ITFAILS;
    if(point_coords[6].size()!=3) ITFAILS;
    if(point_coords[7].size()!=3) ITFAILS;

    if(!soft_equiv(point_coords[0][0],0.0)) ITFAILS;
    if(!soft_equiv(point_coords[0][1],0.0)) ITFAILS;
    if(!soft_equiv(point_coords[0][2],0.0)) ITFAILS;

    if(!soft_equiv(point_coords[1][0],1.0))       ITFAILS;
    if(!soft_equiv(point_coords[1][1],-tan_beta)) ITFAILS;
    if(!soft_equiv(point_coords[1][2],-tan_beta)) ITFAILS;

    if(!soft_equiv(point_coords[2][0],1.0))       ITFAILS;
    if(!soft_equiv(point_coords[2][1],tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[2][2],-tan_beta)) ITFAILS;

    if(!soft_equiv(point_coords[3][0],0.0))       ITFAILS;
    if(!soft_equiv(point_coords[3][1],0.0))       ITFAILS;
    if(!soft_equiv(point_coords[3][2],0.0))       ITFAILS;

    if(!soft_equiv(point_coords[4][0],0.0)) ITFAILS;
    if(!soft_equiv(point_coords[4][1],0.0)) ITFAILS;
    if(!soft_equiv(point_coords[4][2],0.0)) ITFAILS;

    if(!soft_equiv(point_coords[5][0],1.0))       ITFAILS;
    if(!soft_equiv(point_coords[5][1],-tan_beta)) ITFAILS;
    if(!soft_equiv(point_coords[5][2],tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[6][0],1.0))       ITFAILS;
    if(!soft_equiv(point_coords[6][1],tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[6][2],tan_beta)) ITFAILS;

    if(!soft_equiv(point_coords[7][0],0.0))       ITFAILS;
    if(!soft_equiv(point_coords[7][1],0.0))       ITFAILS;
    if(!soft_equiv(point_coords[7][2],0.0))       ITFAILS;

    // check that the mesh returns the correct "next cell"
    if (mesh->next_cell(1,1) !=1) ITFAILS;
    if (mesh->next_cell(1,2) !=0) ITFAILS;
    if (mesh->next_cell(1,3) !=1) ITFAILS;
    if (mesh->next_cell(1,4) !=1) ITFAILS;
    if (mesh->next_cell(1,5) !=1) ITFAILS;
    if (mesh->next_cell(1,6) !=1) ITFAILS;

    // check that mesh can locate a cell give a position
    vector<double> position(3,0.0);
    position[0]=0.5;
    if (mesh->get_cell(position) !=1) ITFAILS;

    position[0]=100;
    if(mesh->get_cell(position)!=-1) ITFAILS;

    // check that the face normals are correct
    vector<double> normal(3,0.0);
    double sin_beta=sin(beta_radians);
    double cos_beta=cos(beta_radians);
    normal=mesh->get_normal(1,1);
    if(!soft_equiv(normal[0],-1.0)) ITFAILS;
    if(!soft_equiv(normal[1], 0.0)) ITFAILS;
    if(!soft_equiv(normal[2], 0.0)) ITFAILS;
    normal=mesh->get_normal(1,2);
    if(!soft_equiv(normal[0], 1.0)) ITFAILS;
    if(!soft_equiv(normal[1], 0.0)) ITFAILS;
    if(!soft_equiv(normal[2], 0.0)) ITFAILS;
    normal=mesh->get_normal(1,3);
    if(!soft_equiv(normal[0], -sin_beta)) ITFAILS;
    if(!soft_equiv(normal[1], -cos_beta)) ITFAILS;
    if(!soft_equiv(normal[2],       0.0)) ITFAILS;
    normal=mesh->get_normal(1,4);
    if(!soft_equiv(normal[0], -sin_beta)) ITFAILS;
    if(!soft_equiv(normal[1],  cos_beta)) ITFAILS;
    if(!soft_equiv(normal[2],       0.0)) ITFAILS;
    normal=mesh->get_normal(1,5);
    if(!soft_equiv(normal[0], -sin_beta)) ITFAILS;
    if(!soft_equiv(normal[1],       0.0)) ITFAILS;
    if(!soft_equiv(normal[2], -cos_beta)) ITFAILS;
    normal=mesh->get_normal(1,6);
    if(!soft_equiv(normal[0], -sin_beta)) ITFAILS;
    if(!soft_equiv(normal[1],       0.0)) ITFAILS;
    if(!soft_equiv(normal[2],  cos_beta)) ITFAILS;

    vector<double> normal_in(3, 0.0);
    for (int face =1; face <=6; face++)
    {
	normal= mesh->get_normal(1,face);
	normal_in=mesh->get_normal_in(1,face);
	for(int dim =0; dim<3;dim++)
	    if(!soft_equiv(normal[dim],-normal_in[dim])) ITFAILS;
    }

    // does the mesh return the correct cell volume
    if (!soft_equiv(mesh->volume(1),          4./3.*tan_beta*tan_beta)) ITFAILS;
    if (!soft_equiv(mesh->get_total_volume(), 4./3.*tan_beta*tan_beta)) ITFAILS;


    // does the mesh return the correct face areas
    if(!soft_equiv(mesh->face_area(1,1),0.0))                  ITFAILS;
    if(!soft_equiv(mesh->face_area(1,2),4.*tan_beta*tan_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(1,3),1./cos_beta))          ITFAILS;
    if(!soft_equiv(mesh->face_area(1,4),1./cos_beta))          ITFAILS;
    if(!soft_equiv(mesh->face_area(1,5),1./cos_beta))          ITFAILS;
    if(!soft_equiv(mesh->face_area(1,6),1./cos_beta))          ITFAILS;

    // check that the mesh returns the correct vertices
    // first, set reference vertices.
    vector<vector<double> > ref_vert;
    ref_vert.resize(3);
    for(int dim =0; dim < 3; dim++)
    {
	ref_vert[dim].resize(8);
	for(int node =0; node<8;node++)
	    ref_vert[dim][node]=0.0;
    }
    //node 1 (lox, loy, loz) - all zeros
    //node 2 (hix, loy, loz)
    ref_vert[0][1]=1.0;
    ref_vert[1][1]=-tan_beta;
    ref_vert[2][1]=-tan_beta;
    //node 3 (hix,hiy,loz)
    ref_vert[0][2]=1.0;
    ref_vert[1][2]=tan_beta;
    ref_vert[2][2]=-tan_beta;
    //node 4 (lox,hiy,loz) - all zeros
    //node 5 (lox, loy, hiz) - all zeros
    //node 6 (hix, loy, hiz)
    ref_vert[0][5]=1.0;
    ref_vert[1][5]=-tan_beta;
    ref_vert[2][5]=tan_beta;
    //node 7 (hix,hiy,hiz)
    ref_vert[0][6]=1.0;
    ref_vert[1][6]=tan_beta;
    ref_vert[2][6]=tan_beta;
    //node 8 (lox,hiy,hiz) - all zeros

    vector<vector<int> > node_num(6,vector<int>(0));
    for (int face =0; face < 6; face++)
	node_num[face].resize(4);
    //face 1 low x
    node_num[0][0]=1;
    node_num[0][1]=4;
    node_num[0][2]=8;
    node_num[0][3]=5;
    //face 2 high x
    node_num[1][0]=2;
    node_num[1][1]=3;
    node_num[1][2]=7;
    node_num[1][3]=6;
    //face 3 low y
    node_num[2][0]=1;
    node_num[2][1]=2;
    node_num[2][2]=6;
    node_num[2][3]=5;
    //face 4 high y
    node_num[3][0]=4;
    node_num[3][1]=3;
    node_num[3][2]=7;
    node_num[3][3]=8;
    //face 5 low z
    node_num[4][0]=1;
    node_num[4][1]=2;
    node_num[4][2]=3;
    node_num[4][3]=4;
    //face 6 high z
    node_num[5][0]=5;
    node_num[5][1]=6;
    node_num[5][2]=7;
    node_num[5][3]=8;

    vector<vector<double> > vertices;
    for(int face=0; face <6; face++)
    {
	vertices = mesh->get_vertices(1,face+1);

	if(vertices.size() !=3)    ITFAILS; // 3 dimensions
	if(vertices[0].size() !=4) ITFAILS; // 4 nodes per face
	if(vertices[1].size() !=4) ITFAILS; // 4 nodes per face
	if(vertices[2].size() !=4) ITFAILS; // 4 nodes per face

	for (int face_node = 0; face_node<4; face_node++)
	    for(int dim =0; dim<3; dim++)
		if(!soft_equiv(vertices[dim][face_node],
			       ref_vert[dim][node_num[face][face_node]-1])) ITFAILS;
    }
    vertices = mesh->get_vertices(1);

    if(vertices.size() !=3)    ITFAILS; // 3 dimensions
    if(vertices[0].size() !=8) ITFAILS; // 8 nodes per cell
    if(vertices[1].size() !=8) ITFAILS; // 8 nodes per cell
    if(vertices[2].size() !=8) ITFAILS; // 8 nodes per cell

    for(int cell_node=0; cell_node <8; cell_node++)
	for(int dim=0; dim <3; dim++)
	    if(!soft_equiv(vertices[dim][cell_node],ref_vert[dim][cell_node])) ITFAILS;

    // checke that the mesh returns the correct face number
    if(mesh->get_bndface("lox",1) != 1) ITFAILS;
    if(mesh->get_bndface("hix",1) != 2) ITFAILS;
    if(mesh->get_bndface("loy",1) != 3) ITFAILS;
    if(mesh->get_bndface("hiy",1) != 4) ITFAILS;
    if(mesh->get_bndface("loz",1) != 5) ITFAILS;
    if(mesh->get_bndface("hiz",1) != 6) ITFAILS;

    // check that, if this were a surface source cell, it would actually be
    // on a vacuum boundary
    vector<int> sur_cell(1,1);
    if(mesh->check_defined_surcells("lox",sur_cell)) ITFAILS;
    if(!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("loy",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("hiy",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("loz",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("hiz",sur_cell)) ITFAILS;

    vector<int> neighbors;
    neighbors = mesh->get_neighbors(1);
    if(neighbors.size() !=6) ITFAILS;
    if(neighbors[0] !=1)     ITFAILS;
    if(neighbors[1] !=0)     ITFAILS;
    if(neighbors[2] !=1)     ITFAILS;
    if(neighbors[3] !=1)     ITFAILS;
    if(neighbors[4] !=1)     ITFAILS;
    if(neighbors[5] !=1)     ITFAILS;

    // make usre the mesh thinks it's not submesh
    if(!mesh->full_Mesh())   ITFAILS;

    // >>>> test the mesh's ability to sample positions. <<<<

    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    if(seed!=1234567)  ITFAILS;
    
    // The explicit-value tests has been hand-checked in a weak sense and,
    // therefore, are more of a regression test.
    vector<double> position_sampled;




    if(rtt_mc_test::passed)
	PASSMSG("Pyramid_Mesh simple_one_cell  problem o.k");
}


//---------------------------------------------------------------------------//
// test the Spyramid_Mesh via the Spyramid_Builder
void build_a_Spyramid()
{
    using rtt_dsxx::SP;
    using rtt_mc_test::Parser;
    using rtt_mc::Spyramid_Builder;
    using rtt_mc::Pyramid_Mesh;
    using std::vector;
    using std::string;
    using rtt_mc::global::pi;
    using std::pow;
    using rtt_dsxx::soft_equiv;
    using std::cos;
    using std::tan;
    using std::sin;
    using std::atan;

    SP<Parser> parser(new Parser("Spyramid_Input"));
    Spyramid_Builder builder(parser);

    // check some of the Spyramid_Builder properties; before build_Mesh, the
    // coordinate system is one-dimensional.  After build_Mesh, it's 3D XYZ

    // test cell region data
    {
	vector<int> regions(5,1);
	regions[2] = 2;
	regions[3] = 2;
	regions[4] = 2;

	if (builder.get_num_regions() !=2)      ITFAILS;
	if (builder.get_regions() != regions)   ITFAILS;
    }

    // test zone mapping
    {
	vector<vector<int> > zone;
	zone.resize(2);
	zone[0].resize(2);
	zone[1].resize(3);

	zone[0][0]=1;
	zone[0][1]=2;
	zone[1][0]=3;
	zone[1][1]=4;
	zone[1][2]=5;

	if(builder.get_num_zones() != 2)           ITFAILS;
	if(builder.get_cells_in_zone(1) !=zone[0]) ITFAILS;
	if(builder.get_cells_in_zone(2) !=zone[1]) ITFAILS;

    }
    //build a Spyramid mesh
    SP<Pyramid_Mesh> mesh= builder.build_Mesh();

    //check defined surface source cells
    {
	vector<vector<int> > ss(1);
	ss[0].resize(1);
	ss[0][0]=5;

	if(builder.get_defined_surcells() !=ss) ITFAILS;

	vector<string> ssp(1);
	ssp[0]="hir";
	if(builder.get_ss_pos() !=ssp) ITFAILS;
    }

    // check zone mapper

    {
	vector<int> zone_field(2);
	zone_field[0] = 1000;
	zone_field[1]= 1001;
	vector<int> cell_field = builder.zone_cell_mapper(zone_field);

	if(cell_field.size() !=mesh->num_cells()) ITFAILS;

	if(cell_field[0] !=1000) ITFAILS;
	if(cell_field[1] !=1000) ITFAILS;
	if(cell_field[2] !=1001) ITFAILS;
	if(cell_field[3] !=1001) ITFAILS;
	if(cell_field[4] !=1001) ITFAILS;
    }

    // check number of cells
    if(builder.num_cells() !=5) ITFAILS;

    // <<<< now do some checks on the mesh itself >>>>
    // (this is a different mesh than the simple_one_cell_Spyramid_Mesh
    

    // check that the mesh returns proper coordinate system information
    if(mesh->get_Coord().get_dim() !=3) ITFAILS;
    if(mesh->get_SPCoord()->get_dim() !=3) ITFAILS;
    if(mesh->get_Coord().get_Coord() !=string("xyz")) ITFAILS;
    if(mesh->get_SPCoord()->get_Coord() != string("xyz")) ITFAILS;

    // check that the mesh returns the proper number of cells
    if (mesh->num_cells() !=5) ITFAILS;

    // calculate r_to_x
    double r_to_x=pow((2.*(1.-cos(pi/4.)))/(tan(pi/4)*tan(pi/4)), 1./3.);

    // check that the x extents are correct for the 1st and last cells
    // the radial boundaries are 0,0.5,1,9/7,13/7,3
    if(!soft_equiv(mesh->get_low_x(1),0.0))           ITFAILS;
    if(!soft_equiv(mesh->get_high_x(1),r_to_x*0.5))   ITFAILS;

    if(!soft_equiv(mesh->get_low_x(5), r_to_x*13./7.)) ITFAILS;
    if(!soft_equiv(mesh->get_high_x(5), r_to_x*3.))    ITFAILS;

    // check that cell midpoints are correct for 1st and last cell
    if(!soft_equiv(mesh->get_x_midpoint(1), 0.5*0.5*r_to_x))  ITFAILS;
    if(!soft_equiv(mesh->get_y_midpoint(1),0.0))              ITFAILS;
    if(!soft_equiv(mesh->get_z_midpoint(1),0.0))              ITFAILS;

    if(!soft_equiv(mesh->get_x_midpoint(5), 0.5*(34./7.)*r_to_x))ITFAILS;
    if(!soft_equiv(mesh->get_y_midpoint(5),0.0))                 ITFAILS;
    if(!soft_equiv(mesh->get_z_midpoint(5),0.0))                 ITFAILS;

    // check that the cell dimensions are correct for 1st and last cells;
    double tan_beta =sqrt(pi)/2.;

    if(!soft_equiv(mesh->dim(1,1),0.5*r_to_x))   ITFAILS;
    if(!soft_equiv(mesh->dim(2,1),0.5*tan_beta*r_to_x)) ITFAILS;
    if(!soft_equiv(mesh->dim(3,1),0.5*tan_beta*r_to_x)) ITFAILS;

    if(!soft_equiv(mesh->dim(1,5),8./7.*r_to_x))  ITFAILS;
    if(!soft_equiv(mesh->dim(2,5),34./7.*tan_beta*r_to_x)) ITFAILS;
    if(!soft_equiv(mesh->dim(3,5),34./7.*tan_beta*r_to_x)) ITFAILS;

    // check that graphics cell type is correct
    vector<int> cell_type=mesh->get_cell_types();
    if(cell_type.size() !=5)                     ITFAILS;
    for(int i=0; i<5; i++)
	if(cell_type[0]!=rtt_viz::eight_node_hexahedron) ITFAILS;

    // check that mesh returns the proper point coordinates (1st & last cell)
    vector<vector<double> > point_coords = mesh->get_point_coord();
    if(point_coords.size() != 40)                       ITFAILS;
    for(int p=0;p<40;p++)
	if(point_coords[p].size() != 3)                ITFAILS;

    // first cell, low z face, counter clockwise from the lower left, facing +z
    if(!soft_equiv(point_coords[0][0], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[0][1], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[0][2], 0.0))           ITFAILS;

    if(!soft_equiv(point_coords[1][0], 0.5*r_to_x))           ITFAILS;
    if(!soft_equiv(point_coords[1][1],-0.5*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[1][2],-0.5*r_to_x*tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[2][0], 0.5*r_to_x))           ITFAILS;
    if(!soft_equiv(point_coords[2][1], 0.5*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[2][2],-0.5*r_to_x*tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[3][0], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[3][1], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[3][2], 0.0))           ITFAILS;

    //first cell, high z face, counter clocwise from the lower left, facing +z
    if(!soft_equiv(point_coords[4][0], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[4][1], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[4][2], 0.0))           ITFAILS;

    if(!soft_equiv(point_coords[5][0], 0.5*r_to_x))           ITFAILS;
    if(!soft_equiv(point_coords[5][1],-0.5*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[5][2], 0.5*r_to_x*tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[6][0], 0.5*r_to_x))           ITFAILS;
    if(!soft_equiv(point_coords[6][1], 0.5*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[6][2], 0.5*r_to_x*tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[7][0], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[7][1], 0.0))           ITFAILS;
    if(!soft_equiv(point_coords[7][2], 0.0))           ITFAILS;


    // last cell, low z face, counter clockwise from the lower left, facing +z
    if(!soft_equiv(point_coords[32][0], 13./7.*r_to_x))           ITFAILS;
    if(!soft_equiv(point_coords[32][1],-13./7.*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(point_coords[32][2],-13./7.*r_to_x*tan_beta))  ITFAILS;

    if(!soft_equiv(point_coords[33][0], 3.*r_to_x))              ITFAILS;
    if(!soft_equiv(point_coords[33][1],-3.*r_to_x*tan_beta))     ITFAILS;
    if(!soft_equiv(point_coords[33][2],-3.*r_to_x*tan_beta))     ITFAILS;
    
    if(!soft_equiv(point_coords[34][0], 3.*r_to_x))              ITFAILS;
    if(!soft_equiv(point_coords[34][1], 3.*r_to_x*tan_beta))     ITFAILS;
    if(!soft_equiv(point_coords[34][2],-3.*r_to_x*tan_beta))     ITFAILS;

    if(!soft_equiv(point_coords[35][0], 13./7.*r_to_x))          ITFAILS;
    if(!soft_equiv(point_coords[35][1], 13./7.*r_to_x*tan_beta)) ITFAILS;
    if(!soft_equiv(point_coords[35][2],-13./7*r_to_x*tan_beta))  ITFAILS;

    // last cell, high z face, counter clockwise from the lower left, facing +z
    if(!soft_equiv(point_coords[36][0], 13./7.*r_to_x))          ITFAILS;
    if(!soft_equiv(point_coords[36][1],-13./7.*r_to_x*tan_beta)) ITFAILS;
    if(!soft_equiv(point_coords[36][2], 13./7.*r_to_x*tan_beta)) ITFAILS;

    if(!soft_equiv(point_coords[37][0], 3.*r_to_x))              ITFAILS;
    if(!soft_equiv(point_coords[37][1],-3.*r_to_x*tan_beta))     ITFAILS;
    if(!soft_equiv(point_coords[37][2], 3.*r_to_x*tan_beta))     ITFAILS;
    
    if(!soft_equiv(point_coords[38][0], 3.*r_to_x))              ITFAILS;
    if(!soft_equiv(point_coords[38][1], 3.*r_to_x*tan_beta))     ITFAILS;
    if(!soft_equiv(point_coords[38][2], 3.*r_to_x*tan_beta))     ITFAILS;

    if(!soft_equiv(point_coords[39][0], 13./7.*r_to_x))          ITFAILS;
    if(!soft_equiv(point_coords[39][1], 13./7.*r_to_x*tan_beta)) ITFAILS;
    if(!soft_equiv(point_coords[39][2], 13./7.*r_to_x*tan_beta)) ITFAILS;

    // check that mesh returns the correct "next cell"
    // low and high y and z should be reflecting
    for(int c=1; c<=5; c++)
    {
	if(mesh->next_cell(c,3)!=c) ITFAILS;
	if(mesh->next_cell(c,4)!=c) ITFAILS;
	if(mesh->next_cell(c,5)!=c) ITFAILS;
	if(mesh->next_cell(c,6)!=c) ITFAILS;
    }
    if(mesh->next_cell(1,1) !=1) ITFAILS;
    if(mesh->next_cell(1,2) !=2) ITFAILS;
    if(mesh->next_cell(2,1) !=1) ITFAILS;
    if(mesh->next_cell(2,2) !=3) ITFAILS;
    if(mesh->next_cell(3,1) !=2) ITFAILS;
    if(mesh->next_cell(3,2) !=4) ITFAILS;
    if(mesh->next_cell(4,1) !=3) ITFAILS;
    if(mesh->next_cell(4,2) !=5) ITFAILS;
    if(mesh->next_cell(5,1) !=4) ITFAILS;
    if(mesh->next_cell(5,2) !=0) ITFAILS;
    
    // chekc that mesh can locate a cell given a position
    vector<double> position(3, 0.0);
    position[0]=0.25*r_to_x;
    if(mesh->get_cell(position) !=1) ITFAILS;

    position[0]= 100.0;
    if(mesh->get_cell(position) !=-1) ITFAILS;

    position[0]=9./7*r_to_x;
    position[1]=r_to_x*tan_beta;
    position[2]=-r_to_x*tan_beta;
    if(mesh->get_cell(position) != 4) ITFAILS;

    // check that the face normals are correct
    vector<double> normal(3,0.0);
    double cos_beta=cos(atan(tan_beta));
    double sin_beta=sin(atan(tan_beta));
    normal=mesh->get_normal(1,1);
    if(!soft_equiv(normal[0],-1.0))  ITFAILS;
    if(!soft_equiv(normal[1], 0.0))  ITFAILS;
    if(!soft_equiv(normal[2], 0.0))  ITFAILS;
    normal=mesh->get_normal(1,2);
    if(!soft_equiv(normal[0], 1.0))  ITFAILS;
    if(!soft_equiv(normal[1], 0.0))  ITFAILS;
    if(!soft_equiv(normal[2], 0.0))  ITFAILS;
    normal=mesh->get_normal(1,3);
    if(!soft_equiv(normal[0],-sin_beta))  ITFAILS;
    if(!soft_equiv(normal[1],-cos_beta))  ITFAILS;
    if(!soft_equiv(normal[2], 0.0))       ITFAILS;
    normal=mesh->get_normal(1,4);
    if(!soft_equiv(normal[0],-sin_beta))  ITFAILS;
    if(!soft_equiv(normal[1], cos_beta))  ITFAILS;
    if(!soft_equiv(normal[2], 0.0))       ITFAILS;
    normal=mesh->get_normal(1,5);
    if(!soft_equiv(normal[0],-sin_beta))  ITFAILS;
    if(!soft_equiv(normal[1], 0.0))       ITFAILS;
    if(!soft_equiv(normal[2],-cos_beta))  ITFAILS;
    normal=mesh->get_normal(1,6);
    if(!soft_equiv(normal[0],-sin_beta))  ITFAILS;
    if(!soft_equiv(normal[1], 0.0))       ITFAILS;
    if(!soft_equiv(normal[2], cos_beta))  ITFAILS;

    // now see that all other cells have the same normals
    for (int c = 1; c<=5; c++)
	for(int f=1; f<=6; f++)
	{
	    normal= mesh->get_normal(1,f);
	    vector<double> check = mesh->get_normal(c,f);
	    for(int d= 0; d<3;d++)
		if(!soft_equiv(normal[d],check[d]))  ITFAILS;
	}

    vector<double> normal_in(3,0.0);
    for(int cell =1; cell <=5; cell++)
	for(int face =1; face<=6; face++)
	{
	    normal = mesh->get_normal(1,face);
	    normal_in = mesh->get_normal_in(1,face);

	    for (int dim = 0; dim< 3; dim++)
		if( !soft_equiv(normal[dim],-normal_in[dim]))  ITFAILS;
	}

    // does the mesh return the correct cell volume?
    if(!soft_equiv(mesh->volume(1), 
		   4./3.*tan_beta*tan_beta*pow(0.5*r_to_x,3)))     ITFAILS;
    if(!soft_equiv(mesh->volume(2),
		   4./3.*tan_beta*tan_beta*
		   (pow(1.*r_to_x,3)-pow(0.5*r_to_x,3))))          ITFAILS;
    if(!soft_equiv(mesh->volume(3),
		   4./3.*tan_beta*tan_beta*
		   (pow(9./7.*r_to_x,3)-pow(1.0*r_to_x,3))))       ITFAILS;
    if(!soft_equiv(mesh->volume(4),
		   4./3.*tan_beta*tan_beta*
		   (pow(13./7.*r_to_x,3)-pow(9./7.*r_to_x,3))))    ITFAILS;
    if(!soft_equiv(mesh->volume(5),
		   4./3.*tan_beta*tan_beta*
		   (pow(3.*r_to_x,3)-pow(13./7.*r_to_x,3))))       ITFAILS;

    // does the mesh return the correct total volume?
    if(!soft_equiv(mesh->get_total_volume(), 
		   4./3.*tan_beta*tan_beta*pow(3.0*r_to_x,3)))     ITFAILS;

    // does the mesh return the correct face area?
    if(!soft_equiv(mesh->face_area(1,1),0.0))                      ITFAILS;
    if(!soft_equiv(mesh->face_area(1,2),
		   4.*pow(0.5*tan_beta*r_to_x,2)))                 ITFAILS;
    if(!soft_equiv(mesh->face_area(1,3),
		   pow(0.5*r_to_x,2)/cos_beta))                    ITFAILS;
    if(!soft_equiv(mesh->face_area(1,4),
		   pow(0.5*r_to_x,2)/cos_beta))                    ITFAILS;
    if(!soft_equiv(mesh->face_area(1,5),
		   pow(0.5*r_to_x,2)/cos_beta))                    ITFAILS;
    if(!soft_equiv(mesh->face_area(1,6),
		   pow(0.5*r_to_x,2)/cos_beta))                    ITFAILS;

    if(!soft_equiv(mesh->face_area(2,1),
		   4.*pow(0.5*tan_beta*r_to_x,2)))                  ITFAILS;
    if(!soft_equiv(mesh->face_area(2,2),
		   4.*pow(1.0*tan_beta*r_to_x,2)))                  ITFAILS;
    if(!soft_equiv(mesh->face_area(2,3),
		   (pow(1.0*r_to_x,2)-pow(0.5*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(2,4),
		   (pow(1.0*r_to_x,2)-pow(0.5*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(2,5),
		   (pow(1.0*r_to_x,2)-pow(0.5*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(2,6),
		   (pow(1.0*r_to_x,2)-pow(0.5*r_to_x,2))/cos_beta)) ITFAILS;

    if(!soft_equiv(mesh->face_area(3,1),
		   4.*pow(1.0*tan_beta*r_to_x,2)))                  ITFAILS;
    if(!soft_equiv(mesh->face_area(3,2),
		   4.*pow(9./7.*tan_beta*r_to_x,2)))                ITFAILS;
    if(!soft_equiv(mesh->face_area(3,3),
		   (pow(9./7.*r_to_x,2)-pow(1.0*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(3,4),
		   (pow(9./7.*r_to_x,2)-pow(1.0*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(3,5),
		   (pow(9./7.*r_to_x,2)-pow(1.0*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(3,6),
		   (pow(9./7.*r_to_x,2)-pow(1.0*r_to_x,2))/cos_beta)) ITFAILS;

    if(!soft_equiv(mesh->face_area(4,1),
		   4.*pow(9./7.*tan_beta*r_to_x,2)))                  ITFAILS;
    if(!soft_equiv(mesh->face_area(4,2),
		   4.*pow(13./7.*tan_beta*r_to_x,2)))                 ITFAILS;
    if(!soft_equiv(mesh->face_area(4,3),
		   (pow(13./7.*r_to_x,2)-pow(9./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(4,4),
		   (pow(13./7.*r_to_x,2)-pow(9./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(4,5),
		   (pow(13./7.*r_to_x,2)-pow(9./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(4,6),
		   (pow(13./7.*r_to_x,2)-pow(9./7.*r_to_x,2))/cos_beta)) ITFAILS;

    if(!soft_equiv(mesh->face_area(5,1),
		   4.*pow(13./7.*tan_beta*r_to_x,2)))                 ITFAILS;
    if(!soft_equiv(mesh->face_area(5,2),
		   4.*pow(3.*tan_beta*r_to_x,2)))                     ITFAILS;
    if(!soft_equiv(mesh->face_area(5,3),
		   (pow(3.*r_to_x,2)-pow(13./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(5,4),
		   (pow(3.*r_to_x,2)-pow(13./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(5,5),
		   (pow(3.*r_to_x,2)-pow(13./7.*r_to_x,2))/cos_beta)) ITFAILS;
    if(!soft_equiv(mesh->face_area(5,6),
		   (pow(3.*r_to_x,2)-pow(13./7.*r_to_x,2))/cos_beta)) ITFAILS;


    // check that the mest returns correctly sized vertices
    vector<vector<double> > vertices;
    for(int cell = 0; cell <5; cell++)
    {
	for(int face =0 ; face<6; face++)
	{
	    vertices = mesh->get_vertices(cell+1,face+1);

	    if(vertices.size() !=3 )   ITFAILS; // 3 dimensions
	    if(vertices[0].size() !=4) ITFAILS; // 4 nodes per face
	    if(vertices[1].size() !=4) ITFAILS; // 4 nodes per face
	    if(vertices[2].size() !=4) ITFAILS; // 4 nodes per face
	}

	vertices = mesh->get_vertices(cell+1);
	if(vertices.size() != 3)   ITFAILS;  // 3 dimensions
	if(vertices[0].size() !=8) ITFAILS;  // 8 nodes per cell
	if(vertices[1].size() !=8) ITFAILS;  // 8 nodes per cell
	if(vertices[2].size() !=8) ITFAILS;  // 8 nodes per cell

    }

    // spot check a few vertices
    vertices = mesh->get_vertices(5,5);
    if(!soft_equiv(vertices[0][0], 13./7.*r_to_x))            ITFAILS;
    if(!soft_equiv(vertices[1][0], -13./7.*r_to_x*tan_beta))  ITFAILS;
    if(!soft_equiv(vertices[2][0], -13./7.*r_to_x*tan_beta))  ITFAILS;
    vertices = mesh->get_vertices(2,2);
    if(!soft_equiv(vertices[0][2], 1.*r_to_x))            ITFAILS;
    if(!soft_equiv(vertices[1][2], 1.*r_to_x*tan_beta))   ITFAILS;
    if(!soft_equiv(vertices[1][2], 1.*r_to_x*tan_beta))   ITFAILS;

    // check that the mesh returns the correct face number
    if(mesh->get_bndface("lox",1) !=1) ITFAILS;
    if(mesh->get_bndface("hix",1) !=2) ITFAILS;
    if(mesh->get_bndface("lor",1) !=1) ITFAILS;
    if(mesh->get_bndface("hir",1) !=2) ITFAILS;
    if(mesh->get_bndface("loy",1) !=3) ITFAILS;
    if(mesh->get_bndface("hiy",1) !=4) ITFAILS;
    if(mesh->get_bndface("loz",1) !=5) ITFAILS;
    if(mesh->get_bndface("hiz",1) !=6) ITFAILS;

    // check that surface source cells would actuall be on a vacuum boundary
    vector<int> sur_cell;
    sur_cell.resize(1);
    sur_cell[0] = 5;
    if(!mesh->check_defined_surcells("hix",sur_cell)) ITFAILS;
    if(!mesh->check_defined_surcells("hir",sur_cell)) ITFAILS;
    sur_cell[0]=1;
    if(mesh->check_defined_surcells("lox",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("lor",sur_cell)) ITFAILS;
    sur_cell.resize(5);
    for(int c =0; c<5 ;c++)
	sur_cell[c] =c+1;
    if(mesh->check_defined_surcells("loy",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("hiy",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("loz",sur_cell)) ITFAILS;
    if(mesh->check_defined_surcells("hiz",sur_cell)) ITFAILS;

    // check that the mesh returns the proper cell list
    vector<int> surface_cell_list;
    surface_cell_list =mesh->get_surcells("hir");
    if(surface_cell_list.size()!=1)   ITFAILS;
    if(surface_cell_list[0] !=5)      ITFAILS;

    // check that the neighbors are correct
    vector<int> neighbors;
    neighbors = mesh->get_neighbors(1);
    if(neighbors.size() != 6) ITFAILS;
    if(neighbors[0] !=1)      ITFAILS;
    if(neighbors[1] !=2)      ITFAILS;
    if(neighbors[2] !=1)      ITFAILS;
    if(neighbors[3] !=1)      ITFAILS;
    if(neighbors[4] !=1)      ITFAILS;
    if(neighbors[5] !=1)      ITFAILS;
    neighbors = mesh->get_neighbors(2);
    if(neighbors.size() != 6) ITFAILS;
    if(neighbors[0] !=1)      ITFAILS;
    if(neighbors[1] !=3)      ITFAILS;
    if(neighbors[2] !=2)      ITFAILS;
    if(neighbors[3] !=2)      ITFAILS;
    if(neighbors[4] !=2)      ITFAILS;
    if(neighbors[5] !=2)      ITFAILS;
    neighbors = mesh->get_neighbors(3);
    if(neighbors.size() != 6) ITFAILS;
    if(neighbors[0] !=2)      ITFAILS;
    if(neighbors[1] !=4)      ITFAILS;
    if(neighbors[2] !=3)      ITFAILS;
    if(neighbors[3] !=3)      ITFAILS;
    if(neighbors[4] !=3)      ITFAILS;
    if(neighbors[5] !=3)      ITFAILS;
    neighbors = mesh->get_neighbors(4);
    if(neighbors.size() != 6) ITFAILS;
    if(neighbors[0] !=3)      ITFAILS;
    if(neighbors[1] !=5)      ITFAILS;
    if(neighbors[2] !=4)      ITFAILS;
    if(neighbors[3] !=4)      ITFAILS;
    if(neighbors[4] !=4)      ITFAILS;
    if(neighbors[5] !=4)      ITFAILS;
    neighbors = mesh->get_neighbors(5);
    if(neighbors.size() != 6) ITFAILS;
    if(neighbors[0] !=4)      ITFAILS;
    if(neighbors[1] !=0)      ITFAILS;
    if(neighbors[2] !=5)      ITFAILS;
    if(neighbors[3] !=5)      ITFAILS;
    if(neighbors[4] !=5)      ITFAILS;
    if(neighbors[5] !=5)      ITFAILS;


    // make sure the mesh thinks it's not a submesh
    if(!mesh->full_Mesh()) ITFAILS;

    // check in cell function
    // cell 3
    vector<double> r(3);
    r[0]=9./7.*r_to_x-1.E-13;
    r[1]=r_to_x*tan_beta;
    r[2]=r_to_x*tan_beta;

    if(!mesh->in_cell(3,r))  ITFAILS;

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
	simple_one_cell_Pyramid();

	// test the Spyramid mesh Builder
	build_a_Spyramid();
    }
    catch (assertion &ass)
    {
	cout << "While testing tstPyramid_Mesh, " << ass.what()
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
	    cout << "**** tstPyramid_Mesh Test: PASSED on " 
		 << node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    cout << "Done testing tstPyramid_Mesh on " << node() << endl;
    
    Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstPyramid_Mesh.cc
//---------------------------------------------------------------------------//
