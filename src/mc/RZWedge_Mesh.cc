//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/RZWedge_Mesh.cc
 * \author Todd J. Urbatsch
 * \date   Wed Apr  5 17:32:24 2000
 * \brief  RZWedge_Mesh implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_RZWedge_Mesh_cc__
#define __mc_RZWedge_Mesh_cc__

#include "RZWedge_Mesh.hh"
#include "Constants.hh"
#include "viz/Ensight_Translator.hh"
#include <iostream>
#include <iomanip>

namespace rtt_mc
{

using std::vector;
using std::endl;
using std::setw;
using std::ios;
using std::vector;
using std::fill;
using std::ostream;
using std::string;
using std::setiosflags;

using global::pi;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief RZWedge_Mesh constructor.

 * The constructor requires complete arguments of the constituent data
 * necessary to build an RZWedge_Mesh.  It is expected that an appropriate
 * builder class will message general data through an interface to build the
 * specific data structures needed by RZWedge_Mesh.

 * \param coord_ coordinate system smart pointer

 * \param layout_ AMR_Layout giving cell connectivity information

 * \param cell_xz_extents_ the xz coordinates of each cell in the following
 * form: [cell][low x]; [cell][hi x]; [cell][lo z]; [cell][hi z]

 * \param theta_degrees_ angle of the wedge in degrees 

 * \param submesh[=false] boolean indicator; true if this is a sub-mesh

 */
RZWedge_Mesh::RZWedge_Mesh(rtt_dsxx::SP<Coord_sys> coord_,
			   AMR_Layout &layout_,
			   vf_double &cell_xz_extents_,
			   double theta_degrees_,
			   bool submesh_)
    : coord(coord_),
      layout(layout_), 
      cell_xz_extents(cell_xz_extents_),
      theta_degrees(theta_degrees_),
      submesh(submesh_)
{
    // Check coordinate system class
    Require (coord);
    Require (coord->get_dim() == 3);
    Require (coord->get_Coord() == std::string("xyz"));

    // make sure that the unfolding angle is positive and not obtuse
    Require ((theta_degrees > 0.0) && (theta_degrees <= 90.0));

    // check that the cell-extents vector has num_cells elements 
    // and weakly check that each cell-element has 4 extent-elements
    Check (cell_xz_extents.size() == layout.num_cells());
    Check (cell_xz_extents[0].size() == 4);
    Check (cell_xz_extents[layout.num_cells()-1].size() == 4);
    
    calc_wedge_angle_data(theta_degrees);
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate wedge angle data once for use throughout calculation
 *
 * \param theta_radians unfolding angle in degrees
 */
void RZWedge_Mesh::calc_wedge_angle_data(const double theta_degrees)
{
    Require ((theta_degrees > 0.0) && (theta_degrees <= 90.0));
    
    theta_radians = theta_degrees * 2.0 * rtt_mc::global::pi / 360.0;
    Check (theta_radians > 0.0);
    Check (theta_radians <= rtt_mc::global::pi/2.0);

    // precalculate heavily used trig functions of unfolding angle
    double half_theta = 0.5 * theta_radians;
    tan_half_theta    = std::tan(half_theta);
    sin_half_theta    = std::sin(half_theta);
    cos_half_theta    = std::cos(half_theta);

}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR IMC
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate minimun distance to boundary in an RZWedge_Mesh cell
 * 
 * \param r position
 * \param omega direction
 * \param cell cell containing position r
 * 
 * \return get_db minimun distance to boundary
 * \return face (in arguments) face corresponding to min dist-to-bndry
 */
double RZWedge_Mesh::get_db(const sf_double &r, const sf_double &omega, 
			    int cell, int &face) const
{
    using global::dot;

    Require (r.size() == 3);
    Require (omega.size() == 3);
    Require (global::soft_equiv(dot(omega,omega), 1.0, 1.0e-8));

    // set up 6 dists-to-bndry, initialize to huge value
    // -- always 6 faces in an RZWedge mesh cell. 
    vector<double> distance(6, global::huge);

    // low x face (inner radial face) (internal index 0)
    if (omega[0] < 0.0)
	distance[0] = (get_low_x(cell) - r[0]) / omega[0];

    // high x face (outer radial face) (internal index 1)
    if (omega[0] > 0.0)
	distance[1] = (get_high_x(cell) - r[0]) / omega[0];

    // low y face (-theta/2) (internal index 2) (tan(theta/2)x + y = 0)
    if (dot(omega,get_normal(cell,3)) > 0.0)
	distance[2] = -(r[1] + r[0]*tan_half_theta) /
	    (omega[1] + omega[0]*tan_half_theta);

    // high y face (theta/2) (internal index 3) (tan(theta/2)x - y = 0)
    if (dot(omega,get_normal(cell,4)) > 0.0)
	distance[3] = (r[1] - r[0]*tan_half_theta) /
	    (-omega[1] + omega[0]*tan_half_theta);

    // low z face - (internal index 4)
    if (omega[2] < 0.0)
	distance[4] = (get_low_z(cell) - r[2]) / omega[2];

    // high z face - (internal index 5)
    if (omega[2] > 0.0)
	distance[5] = (get_high_z(cell) - r[2]) / omega[2];

    // find face index and value of minimum(distance[face])
    double min_dist = global::huge;
    int face_index = 0;
    for (int f = 0; f < 6; f++)
	if (distance[f] < min_dist)
	{
	    min_dist   = distance[f];
	    face_index = f;
	}

    // set face index (external: in [1,6])
    Ensure (face_index >= 0 && face_index < 6);
    face = face_index + 1;

    // check that return quantities are within expected limits
    Ensure (min_dist >= 0.0);
    Ensure (face > 0 && face <= 6);
    
    // return minimum distance to boundary 
    return min_dist;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find cell in RZWedge_Mesh given position
 * 
 * \param r position
 * 
 * \return cell cell containing position r
 */
int RZWedge_Mesh::get_cell(const sf_double &r) const
{
    // This is a butt-slow algorithm to merely provide the capability of
    // locating a cell by position.  We do this inefficeint algorithm
    // to allow for an as of yet undetermined numbering for an AMR mesh.
    // We can build in logic later to take advantage of particular
    // number schemes.

    Require (!submesh);

    // initialize cell, flag indicating that cell was located
    int located_cell = 0;
    bool found       = false;

    // find cell that the contains the x- and z-positions
    while (!found)
    {
	located_cell += 1;

	double lox = get_low_x(located_cell);
	double hix = get_high_x(located_cell);
	double loz = get_low_z(located_cell);
	double hiz = get_high_z(located_cell);

	if (r[0] >= lox && r[0] <= hix)
	    if (r[2] >= loz && r[2] <= hiz)
		found = true;

	if (located_cell == num_cells() && !found)
	{
	    located_cell = -1;
	    found = true;
	}
    } 

    // check that the y-position is within the wedge
    if (located_cell > 0)
        Check (std::fabs(r[1]) <= r[0] * tan_half_theta);

    return located_cell;
}

//---------------------------------------------------------------------------//
// return the face number for a given cell boundary (independent of cell)

int RZWedge_Mesh::get_bndface(std_string boundary, int cell) const
{
    // return value
    int face;

    if (boundary == "lox")
	face = 1;
    else if (boundary == "hix")
	face = 2;
    else if (boundary == "loy")
	face = 3;
    else if (boundary == "hiy")
	face = 4;
    else if (boundary == "loz")
	face = 5;
    else if (boundary == "hiz")
	face = 6;
    else
	Insist (0, "Unknown boundary string used!");

    // return the face number given the cell boundary
    return face;
}

//---------------------------------------------------------------------------//
// return a list of cells along a specified boundary

RZWedge_Mesh::sf_int RZWedge_Mesh::get_surcells(std::string boundary) const
{
    Require (!submesh);
    Require (coord->get_dim() == 3);

    // make return vector containing a list of cells along specified boundary
    vector<int> return_list;

    // verify our assumption that cell 1 is in the low x, low z position
    Insist ((layout(1,1).size() == 1) && (layout(1,5).size() == 1) &&
	    (layout(1,3).size() == 1) && (layout(1,4).size() == 1),
	    "Cell 1 should not see refinement on outer boundary and sides!");  
    Insist ((layout(1,1,1)==1) && (layout(1,3,1)==1) && (layout(1,4,1)==1), 
	    "Cell 1 is not reflecting on inside and sides!");
    Insist ((layout(1,5,1) == 0) || (layout(1,5,1) == 1), 
	    "Cell 1 does not have an appropriate outer b.c.!");

    // calculate the cells along the...
    // ...the high r boundary
    if ((boundary == "hix") || (boundary == "hir"))
    {
	// first, work our way along the low z edge to the (high x, low z) cell 
	int start_cell = 1;
	while (layout(start_cell,2,1) != start_cell && 
	       layout(start_cell,2,1) != 0)
	{
	    // get the next cell along the low z edge
	    int next_cell = layout(start_cell,2,1);

	    // ridiculous, unnecessary checks on the next cell 
	    Check (next_cell != 0 && next_cell != start_cell);
	    Check (next_cell > 0 && next_cell <= layout.num_cells());

	    // check that the cell has reflecting wedge attributes
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell)); 

	    // check that the low z edge is on a boundary
	    Check ((layout(next_cell,5,1) == next_cell) ||
		   (layout(next_cell,5,1) == 0));

	    // the next cell is okay
	    start_cell = next_cell;
	}

	// we have the high x starting cell; start filling surface cell list
	int current_cell = start_cell;
	return_list.push_back(current_cell);

	// how many neighbors are there to the high z side?
	int num_cells_hiz = layout.num_cells_across(current_cell,6);
	Check (num_cells_hiz == 1 || num_cells_hiz == 2);

	// keep including the next cell to the high x direction until b.c.
	while (layout(current_cell,2,num_cells_hiz) != current_cell  &&
	       layout(current_cell,2,1) != 0) 
	{
	    // next cell to put into surface cell list
	    int next_cell = layout(current_cell,6,num_cells_hiz);

	    // check that the cell is not outside the system (moot)
	    Check (next_cell != 0 && next_cell != current_cell);

	    // is the next cell a reasonable cell number?
	    Check (next_cell > 0 && next_cell <= layout.num_cells());

	    // check that the next cell is on the high r boundary
	    Check ((layout(next_cell,2,1) == next_cell) ||
		   (layout(next_cell,2,1) == 0)); 

	    // check that the cell has reflecting wedge attributes
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell)); 

	    // accept the next cell and put it into the list
	    current_cell = next_cell;
	    return_list.push_back(current_cell);

	    // get the number of neighbors on the high z side
	    num_cells_hiz = layout.num_cells_across(current_cell,6);
	    Check ((num_cells_hiz == 1) || (num_cells_hiz == 2));
	}

	// check the size of the surface cell list
	Check (return_list.size() > 0);
	Check (return_list.size() <= layout.num_cells()); 
    }

    // ...the low z boundary
    else if (boundary == "loz")
    {
	// include cell 1
	int current_cell = 1;
	return_list.push_back(current_cell);

	// keep including the next cell to the high x direction until b.c.
	while (layout(current_cell,2,1) != current_cell &&
	       layout(current_cell,2,1) != 0) 
	{
	    // next cell to put into surface cell list
	    int next_cell = layout(current_cell,2,1);

	    // check that the cell is not outside the system (moot)
	    Check (next_cell != 0 && next_cell != current_cell);

	    // is the next cell a reasonable cell number?
	    Check (next_cell > 0 && next_cell <= layout.num_cells());

	    // check that the next cell is on the low z boundary
	    Check ((layout(next_cell,5,1) == next_cell) ||
		   (layout(next_cell,5,1) == 0)); 

	    // check that the cell has reflecting wedge attributes
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell)); 

	    // accept the next cell and put it into the list
	    current_cell = next_cell;
	    return_list.push_back(current_cell);
	}

	// check the size of the surface cell list
	Check (return_list.size() > 0);
	Check (return_list.size() <= layout.num_cells()); 
    }

    // ...the high z boundary
    else if (boundary == "hiz")
    {
	// first, work our way along the low x edge to the low x, high z cell 
	int start_cell = 1;
	while (layout(start_cell,6,1) != start_cell && 
	       layout(start_cell,6,1) != 0)
	{
	    // get the next cell along the low x edge
	    int next_cell = layout(start_cell,6,1);

	    // ridiculous, unnecessary checks on the next cell 
	    Check (next_cell != 0 && next_cell != start_cell);
	    Check (next_cell > 0 && next_cell <= layout.num_cells());

	    // check that the cell has reflecting wedge attributes
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell) &&
		   (layout(next_cell,1,1) == next_cell)); 

	    // the next cell is okay
	    start_cell = next_cell;
	}

	// we have the high z starting cell; start filling surface cell list
	int current_cell = start_cell;
	return_list.push_back(current_cell);

	// how many neighbors are there to the high x side?
	int num_cells_hir = layout.num_cells_across(current_cell,2);
	Check (num_cells_hir == 1 || num_cells_hir == 2);

	// keep including the next cell to the high x direction until b.c.
	while (layout(current_cell,2,num_cells_hir) != current_cell  &&
	       layout(current_cell,2,1) != 0) 
	{
	    // next cell to put into surface cell list
	    int next_cell = layout(current_cell,2,num_cells_hir);

	    // check that the cell is not outside the system (moot)
	    Check (next_cell != 0 && next_cell != current_cell);

	    // is the next cell a reasonable cell number?
	    Check (next_cell > 0 && next_cell <= layout.num_cells());

	    // check that the next cell is on the high z boundary
	    Check ((layout(next_cell,6,1) == next_cell) ||
		   (layout(next_cell,6,1) == 0)); 

	    // check that the cell has reflecting wedge attributes
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell)); 

	    // accept the next cell and put it into the list
	    current_cell = next_cell;
	    return_list.push_back(current_cell);

	    // get the number of neighbors on the high x side
	    num_cells_hir = layout.num_cells_across(current_cell,2);
	    Check ((num_cells_hir == 1) || (num_cells_hir == 2));
	}

	// check the size of the surface cell list
	Check (return_list.size() > 0);
	Check (return_list.size() <= layout.num_cells()); 
	    
    }
    else
	Insist(0, "Unknown or invalid (lor/x,loy,hiy) surf in get_surcells!") 

    // return vector
    return return_list;
}

//---------------------------------------------------------------------------//
// check that a user-/host-defined set of surface source cells actually
// resides on the surface of the system (requires a vacuum bnd).

bool RZWedge_Mesh::check_defined_surcells(const std_string ss_face,
					  const sf_int &ss_list) const
{
    // a weak check on number of surface cells
    Check (ss_list.size() <= num_cells());

    for (int ss_indx = 0; ss_indx < ss_list.size(); ss_indx++)
    {
        // convert face on which ss resides from string to int.
        // despite its args, get_bndface actually has no cell dependence
        int ss_face_num = get_bndface(ss_face, ss_list[ss_indx]);

        // get bnd condition on ss face; had better be vacuum (0)
        int bc = layout(ss_list[ss_indx], ss_face_num, 1);
        if (bc != 0)
            return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the vertices for a given RZWedge_Mesh cell. 
 *
 * During normal use of the RZWedge_Mesh, the explicit cell vertices are not
 * required.  However, they are needed for graphics dumps.  The RZWedge_Mesh
 * is always in 3-D XYZ geometry and always has eight vertices (each cell on the
 * radial axis has two sets of conincident vertices).
 * 
 * \param cell RZWedge_Mesh (global) cell number
 * 
 * \return vertices - the coordinates of the cell's eight vertices 
 */
RZWedge_Mesh::vf_double RZWedge_Mesh::get_vertices(int cell) const
{
    Require (cell > 0 && cell <= num_cells());
    Require (coord->get_dim() == 3);

    const int num_verts_face = 4;
    const int num_verts_cell = 8;
    const int loz_face       = 5;
    const int hiz_face       = 6;
	
    // get the vertices for the low z face of the cell
    vf_double cell_vertices = get_vertices(cell, loz_face);

    // get the vertices for the high z face of the cell
    vf_double hiz_face_vertices = get_vertices(cell, hiz_face);
    
    Require (cell_vertices.size() == coord->get_dim());
    Require (cell_vertices[0].size() == num_verts_face);
    Require (cell_vertices[1].size() == num_verts_face);
    Require (cell_vertices[2].size() == num_verts_face);

    Require (hiz_face_vertices.size() == coord->get_dim());
    Require (hiz_face_vertices[0].size() == num_verts_face);
    Require (hiz_face_vertices[1].size() == num_verts_face);
    Require (hiz_face_vertices[2].size() == num_verts_face);

    // add the hiz face vertices to loz face vertices
    for (int v = 0; v < coord->get_dim(); v++)
	cell_vertices[v].insert(cell_vertices[v].end(),
				hiz_face_vertices[v].begin(),
				hiz_face_vertices[v].end()); 

    // make sure each dimension has 8 entries -- one per vertex
    Ensure (cell_vertices[0].size() == num_verts_cell);
    Ensure (cell_vertices[1].size() == num_verts_cell);
    Ensure (cell_vertices[2].size() == num_verts_cell);

    // return the cell vertices
    return cell_vertices;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the vertices for a given face of an RZWedge_Mesh cell. 
 *
 * During normal use of the RZWedge_Mesh, the explicit cell vertices are not
 * required.  However, they are needed for graphics dumps.  The RZWedge_Mesh
 * is always in 3-D XYZ geometry and each face always has four vertices (each
 * cell on the radial axis has two sets of coincident vertices).
 * 
 * \param cell RZWedge_Mesh (global) cell number
 * \param face face of cell
 * 
 * \return vertices the coordinates of the four vertices defining a face
 */
RZWedge_Mesh::vf_double RZWedge_Mesh::get_vertices(int cell, int face) const
{
    Require (face > 0 && face <= 6);
    Require (cell > 0 && cell <= layout.num_cells());
    Require (coord->get_dim() == 3);

    vf_double face_vertices(coord->get_dim());
    sf_double single_vert(coord->get_dim(), 0.0);
    
    double lox     = get_low_x(cell);
    double hix     = get_high_x(cell);
    double small_y = lox * tan_half_theta; 
    double large_y = hix * tan_half_theta;
    double loz     = get_low_z(cell);
    double hiz     = get_high_z(cell);
    
    // low x face or high x face
    if (face == 1 || face == 2)
    {
	double y_variable;
	if (face == 1)
	{
	    single_vert[0] = lox;
	    y_variable = small_y;
	}
	else if (face == 2)
	{
	    single_vert[0] = hix;
	    y_variable = large_y;
	}

	single_vert[1] = -y_variable;
	single_vert[2] = loz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[1] = y_variable;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[2] = hiz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[1] = -y_variable;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);
    }

    // low y face or high y face
    else if (face == 3 || face ==4)
    {
	double y_plusminus;
	if (face == 3)
	    y_plusminus = -1.0;
	else if (face == 4)
	    y_plusminus = 1.0;

	single_vert[0] = lox;
	single_vert[1] = y_plusminus * small_y;
	single_vert[2] = loz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] = hix;
	single_vert[1] = y_plusminus * large_y;
	single_vert[2] = loz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] = hix;
	single_vert[1] = y_plusminus * large_y;
	single_vert[2] = hiz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] = lox;
	single_vert[1] = y_plusminus * small_y;
	single_vert[2] = hiz;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);
    }

    // low z face  or  high z face
    else if (face == 5 || face == 6)
    {
	if (face == 5)
	    single_vert[2] = loz;
	else if (face == 6)
	    single_vert[2] = hiz;

	single_vert[0] = lox;
	single_vert[1] = -small_y;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);
	    
	single_vert[0] = hix;
	single_vert[1] = -large_y;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] = hix;
	single_vert[1] = large_y;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] = lox;
	single_vert[1] = small_y;
	for (int i = 0; i < coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);
    }

    // check that there are 3 dimensions and 4 vertices
    Ensure (face_vertices.size() == coord->get_dim());
    Ensure (face_vertices[0].size() == 4);
    Ensure (face_vertices[1].size() == 4);
    Ensure (face_vertices[2].size() == 4);
   
    // return the four vertices for the face
    return face_vertices;
}

//---------------------------------------------------------------------------//
// functions required for graphics dumps
//---------------------------------------------------------------------------//
/*!
 * \brief Return the cell type for each cell in the RZWedge_Mesh
 */
RZWedge_Mesh::sf_int RZWedge_Mesh::get_cell_types() const
{
    vector<int> cell_type(layout.num_cells());

    // all cells in an RZWedge_Mesh are general, 8-node hexedrons
    std::fill(cell_type.begin(), cell_type.end(),
	      rtt_viz::eight_node_hexahedron); 

    return cell_type;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Return the coordinates of all nodes in the mesh.
 *
 * For each cell in the RZWedge_Mesh, the coordinates of all eight nodes are
 * returned.  Thus, in an unrefined RZWedge_Mesh, all interior nodes are
 * replicated eight times.  This approach, although memory-inefficient, is
 * easier to code especially considering that the RZWedge_Mesh does not 
 * already have the raw vertex data.
 *
 */
RZWedge_Mesh::vf_double RZWedge_Mesh::get_point_coord() const
{
    // number of vertices per cell is always 8; always 3D
    const int num_verts_cell = 8;
    double vert_index;
    Check (coord->get_dim() == 3);

    // weakly check the validity of num_cells()
    Check (num_cells() > 0);

    // initialize the return vector  
    vector<vector<double> > return_coord(num_cells() * num_verts_cell);

    // for each cell, get vertices and reverse the columns and rows
    for (int cell = 1; cell <= num_cells(); cell++)
    {
	// get the cell's vertices[dim][8]
	vector<vector<double> > cell_verts =
	    RZWedge_Mesh::get_vertices(cell); 

	// check validity of cell vertices vector
	Check (cell_verts.size() == coord->get_dim());
	
	// loop over all 8 nodes for this cell
	for (int node = 0; node < num_verts_cell; node++)
	{
	    // calculate running vertex index
	    vert_index = (cell-1)*num_verts_cell + node;

	    // resize each vertice's return_coord to num dimensions
	    return_coord[vert_index].resize(coord->get_dim());

	    // re-assign point coordinates to return vector
	    for (int dim = 0; dim < coord->get_dim(); dim++)
	    {
		Check (cell_verts[dim].size() == num_verts_cell);
		return_coord[vert_index][dim] = cell_verts[dim][node];
	    }
	}
    }

    // return the coordinates for all nodes of all cells
    return return_coord;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get cell-pairing data that matches the coordinate data returned by
 * point coord.
 */
RZWedge_Mesh::vf_int RZWedge_Mesh::get_cell_pair() const
{
    // each cell points to eight vertices, this is ALWAYS a 3D mesh
    vf_int cp(layout.num_cells(), sf_int(8));

    // the coordinates are cyclic starting from the low z to the high z face
    // going clockwise in the xy plane
    int counter = 0;
    for (int cell = 0; cell < cp.size(); cell++)
    {
	for (int node = 0; node < cp[cell].size(); node++)
	    cp[cell][node] = ++counter;

	Check (((cell+1) * counter) % 8 == 0);
    }

    Ensure (counter == cp.size() * 8);
    return cp;
}

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//
// print out the whole mesh

void RZWedge_Mesh::print(ostream &out) const
{
    out << endl;
    out << ">>> MESH <<<" << endl;
    out << "============" << endl;
    out << endl;
    
    out.precision(4);
    out << "--- Wedge Angle --" << endl;
    out << setw(10) << setiosflags(ios::fixed) 
	<< theta_degrees << " degrees" << endl;
    out << setw(10) << setiosflags(ios::fixed)
	<< theta_radians << " radians" << endl;
    out << "------------------" << endl;
    out << endl;

    for (int cell = 1; cell <= num_cells(); cell++)
	print(out, cell);
}

//---------------------------------------------------------------------------//
// print individual cells

void RZWedge_Mesh::print(ostream &output, int cell) const
{
    // print out content info for one cell
    output << "+++++++++++++++" << endl;
    output << "---------------" << endl;
    output << "Cell : "         << cell << endl;
    output << "---------------" << endl;
    output << "Dimensions "     << endl;
    output << "---------------" << endl;

    output.precision(4);
    output << setw(10) << setiosflags(ios::scientific) 
	   << " x range:  [" << get_low_x(cell) << " , " 
	   << get_high_x(cell) << "];  dx : " 
	   << get_high_x(cell)-get_low_x(cell) << endl;
    output << setw(10) << setiosflags(ios::scientific) 
	   << " z range:  [" << get_low_z(cell) << " , " 
	   << get_high_z(cell) << "];  dz : " 
	   << get_high_z(cell)-get_low_z(cell) << endl;

    output << "---------------" << endl;
    output << "Layout "         << endl;
    output << "---------------" << endl;
    layout.print(output, cell);
    output << "+++++++++++++++" << endl;
}

//---------------------------------------------------------------------------//
// Overloaded operators
//---------------------------------------------------------------------------//
// overloaded == for design-by-contract

bool RZWedge_Mesh::operator==(const RZWedge_Mesh &rhs) const
{
    using rtt_mc::global::soft_equiv;

    // check to see that the Layouts are equal
    if (layout != rhs.layout)
        return false;

    // check the XZ extents of the cells
    for (int cell = 1; cell <= layout.num_cells(); cell++)
    {
	if (!soft_equiv(get_low_x(cell),rhs.get_low_x(cell))) 
	    return false;
	if (!soft_equiv(get_high_x(cell),rhs.get_high_x(cell))) 
	    return false;
	if (!soft_equiv(get_low_z(cell),rhs.get_low_z(cell))) 
	    return false;
	if (!soft_equiv(get_high_z(cell),rhs.get_high_z(cell))) 
	    return false;
    }	    

    // if we haven't returned, then the two meshes must be equal
    return true;
}

//---------------------------------------------------------------------------//
// overloaded output operator

std::ostream& operator<<(std::ostream &output, const RZWedge_Mesh &object)
{
    object.print(output);
    return output;
}


} // end namespace rtt_mc

#endif                          // __mc_RZWedge_Mesh_cc__

//---------------------------------------------------------------------------//
//                        end of mc/RZWedge_Mesh.cc
//---------------------------------------------------------------------------//
