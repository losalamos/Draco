//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/RZWedge_Mesh.cc
 * \author Todd J. Urbatsch
 * \date   Wed Apr  5 17:32:24 2000
 * \brief  RZWedge_Mesh implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_RZWedge_Mesh_cc
#define rtt_mc_RZWedge_Mesh_cc

#include "RZWedge_Mesh.hh"
#include "XYZCoord_sys.hh"
#include "Constants.hh"
#include "Math.hh"
#include "viz/Ensight_Translator.hh"
#include "ds++/Packing_Utils.hh"
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <algorithm>

namespace rtt_mc
{

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

    // precalculate heavily used angle quantities    
    calc_wedge_angle_data(theta_degrees);

    // Calculate and assign the on-processor total volume
    calc_total_volume();
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
/*!
 * \brief Calculate and set the total (on-processor) volume of the mesh.
 *
 * \return total volume of on-processor RZWedge cells
 */
void RZWedge_Mesh::calc_total_volume()
{
    Require (num_cells() > 0);

    // initialize private data member
    total_volume = 0.0;

    // sum local cell volumes
    for (int cell = 1; cell <= num_cells(); cell++)
	total_volume += volume(cell);

    Ensure (total_volume > 0.0);
}

//---------------------------------------------------------------------------//
// INTERFACE NOT SPECIFIC FOR IMC
//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a position vector is within a cell.

 * \param cell cell index

 * \param r sf_double of position (3 dim vector)

 * \return true if r is in cell; false otherwise

 */
bool RZWedge_Mesh::in_cell(int cell, const sf_double &r) const
{
    using rtt_mc::global::soft_equiv;

    Require (r.size() == 3);
    Require (cell > 0 && cell <= layout.num_cells());

    // first check x dimension
    if ((r[0] < cell_xz_extents[cell-1][0] &&
	 !soft_equiv(r[0], cell_xz_extents[cell-1][0]))
	|| 
	(r[0] > cell_xz_extents[cell-1][1] &&
	 !soft_equiv(r[0], cell_xz_extents[cell-1][1])))
	return false;

    // check z dimension
    if ((r[2] < cell_xz_extents[cell-1][2] && 
	 !soft_equiv(r[2], cell_xz_extents[cell-1][2]))
	||
	(r[2] > cell_xz_extents[cell-1][3] &&
	 !soft_equiv(r[2], cell_xz_extents[cell-1][3])))
	return false;

    // check y dimension (x cannot be negative)
    if ((r[1] < -(r[0] * tan_half_theta) &&
	 !soft_equiv(r[1], -(r[0] * tan_half_theta)))
	||
	(r[1] > (r[0] * tan_half_theta) &&
	 !soft_equiv(r[1], (r[0] * tan_half_theta))))
	return false;
    
    return true;
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
    using std::vector;
    using global::dot;

    Check (r.size() == 3);
    Check (omega.size() == 3);
    Check (global::soft_equiv(dot(omega,omega), 1.0, 1.0e-5));

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

    if (boundary == "lox" || boundary == "lor")
	face = 1;
    else if (boundary == "hix" || boundary == "hir")
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
    using std::vector;

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
	    Check (next_cell > 0  && next_cell <= layout.num_cells());

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
	while (layout(current_cell,6,num_cells_hiz) != current_cell  &&
	       layout(current_cell,6,1) != 0) 
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
/*! 
 * \brief Sample a position on the surface of a sphere inside a cell.
 *
 * This function is required by the IMC_MT concept.
 *
 * As opposed to the OS_Mesh implementation of this function, this function
 * actually tracks a distance equal to radius.  The normal is the direction
 * of the ray when reaching that distance.
 * 
 * \param cell cell index
 * \param origin sphere origin
 * \param radius sphere radius
 * \param random random number object
 * \return a pair of vector<double> where the first element of the pair is
 * the position on the surface of the sphere and the second element of the
 * pair is the normal of the sphere at that position
 */
RZWedge_Mesh::pair_sf_double RZWedge_Mesh::sample_pos_on_sphere(
    int              cell, 
    const sf_double &origin,
    double           radius,
    rng_Sprng       &random) const
{
    Require (cell > 0);
    Require (cell <= layout.num_cells());
    Require (origin.size() == 3);
    Require (in_cell(cell, origin));

    // checks to make sure sphere is in cell in x and z dimensions
    Require (origin[0] - radius >= get_low_x(cell));
    Require (origin[0] + radius <= get_high_x(cell));
    Require (origin[2] - radius >= get_low_z(cell));
    Require (origin[2] + radius <= get_high_z(cell));

    // get initial position and direction to track
    sf_double r     = origin;
    sf_double omega = coord->sample_isotropic_dir(random);

    // distance we have to track
    double track = radius;
    
    // track until we have gone the radial distance
    double d_bnd     = 0.0;
    int face         = 0;
    double factor    = 0.0;
    sf_double normal;
    while (track > 0.0)
    {
	// determine shortest distance to boundary
	d_bnd = get_db(r, omega, cell, face);
	Check (face == 3 || face == 4);

	// process a reflection on a y face
	if (d_bnd < track)
	{
	    // stream to the face
	    r[0] = r[0] + d_bnd * omega[0];
	    r[1] = r[1] + d_bnd * omega[1];
	    r[2] = r[2] + d_bnd * omega[2];

	    // adjust the remaining track length
	    track -= d_bnd;
	    Check (track >= 0.0);

	    // calculate face normal
	    normal = get_normal(cell, face);
	    Check (normal.size() == 3);
	    
	    // specularly reflect angle
	    factor    = rtt_mc::global::dot(omega, normal);
	    omega[0] -= 2.0 * factor * normal[0];
	    omega[1] -= 2.0 * factor * normal[1];
	    omega[2] -= 2.0 * factor * normal[2];
	}
	
	// stream until we are finished
	else
	{
	    r[0] = r[0] + track * omega[0];
	    r[1] = r[1] + track * omega[1];
	    r[2] = r[2] + track * omega[2];

	    // we are finished tracking
	    track = 0.0;
	}
    }

    // assign and return position and normal, the normal is the last
    // direction the ray had before reaching the specfied track distance 
    pair_sf_double pos_and_norm = std::make_pair(r, omega);

    // checks; we do not check that omega is properly normalized here because
    // it is checked in the preceding and ensuing transport steps
    Ensure (in_cell(cell, pos_and_norm.first));

    // return
    return pos_and_norm;
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
// Interface for graphics dumps
//---------------------------------------------------------------------------//
/*!
 * \brief Return the cell type for each cell in the RZWedge_Mesh
 */
RZWedge_Mesh::sf_int RZWedge_Mesh::get_cell_types() const
{
    std::vector<int> cell_type(layout.num_cells());

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
    using std::vector;

    // number of vertices per cell is always 8; always 3D
    const int num_verts_cell = 8;
    int vert_index;
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

void RZWedge_Mesh::print(std::ostream &out) const
{
    using std::setw;
    using std::setiosflags;
    using std::ios;
    using std::endl;

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

void RZWedge_Mesh::print(std::ostream &output, int cell) const
{
    using std::endl;
    using std::setiosflags;
    using std::setw;
    using std::ios;

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

//---------------------------------------------------------------------------//
// Mesh Packing Interface
//---------------------------------------------------------------------------//
/*!

 * \brief Pack up a mesh into a Pack struct for communication and
 * persistence.

 * The cell list provides the cells to pack up.  It also is a map from a full
 * mesh to a spatially decomposed mesh.  Thus, the packed mesh will only
 * contain cells in the cell list with the provided mappings.

 * The packer will not produce an "exact" copy of the mesh even if the
 * current_mesh_to_new_mesh mapping is one to one.  It will produce an
 * equivalent copy (the internal data will be organized differently), and
 * operator== will fail on such a comparison.  To produce an exact copy, call
 * pack without any arguments.
 
 * \param current_mesh_to_new_mesh list of cells to include in the packed
 * mesh, set this to NULL to produce an exact copy

 */
RZWedge_Mesh::SP_Pack RZWedge_Mesh::pack(const sf_int &current_to_new) const
{
    Require (current_to_new.size() == layout.num_cells() ||
	     current_to_new.size() == 0);

    // determine whether this is exact replication or sub-packing
    sf_int current_to_new_replicate;
    bool   replicate;
    if (current_to_new.size() == 0)
    {
	replicate = true;
	current_to_new_replicate.resize(layout.num_cells());
	
	for (int cell = 1; cell <= layout.num_cells(); cell++)
	    current_to_new_replicate[cell-1] = cell;
    }
    else 
    {
	replicate = false;
    }
    
    // the coordinate system is always XYZ
    Ensure (typeid(*coord) == typeid(XYZCoord_sys));

    // packup layout data
    AMR_Layout::SP_Pack packed_layout;
    int                 layout_size;
    int                 num_packed_cells;

    // packed extents data
    int   extent_size = 0;
    char *extent_data = 0;

    // pack up the mesh
    if (replicate)
    {
	// packup the layout
	packed_layout = layout.pack(current_to_new_replicate);
	layout_size   = packed_layout->get_size();
	Check (layout_size >= 1);

	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells == layout.num_cells());

	// calculate extent size
	extent_size = (4 * num_packed_cells) * sizeof(double);
	extent_data = new char[extent_size];

	// pack extents
	pack_extents(current_to_new_replicate, extent_data, extent_size,
		     num_packed_cells);
    }
    else
    {
	// packup the layout
	packed_layout = layout.pack(current_to_new);
	layout_size   = packed_layout->get_size();
	Check (layout_size >= 1);

	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells <= layout.num_cells());

	// calculate extent size
	extent_size = (4 * num_packed_cells) * sizeof(double);
	extent_data = new char[extent_size];
	
	// pack extents
	pack_extents(current_to_new, extent_data, extent_size,
		     num_packed_cells);
    }
    Check (num_packed_cells >= 0 && num_packed_cells <= layout.num_cells());
    Check (extent_size == (num_packed_cells * 4) * sizeof(double));
    Check (extent_data != 0);

    // now pack up the mesh
    
    // ints (1-size of packed layout; 1-num_packed_cells; packed layout)
    int total_ints    = (2 + layout_size) * sizeof(int);
    
    // doubles (1-theta angle)
    int total_doubles = 1 * sizeof(double);

    // chars
    int total_chars   = extent_size;

    // allocate space
    int   size = total_ints + total_doubles + total_chars;
    char *data = new char[size];

    // pack up the mesh
    
    // make a packer object
    rtt_dsxx::Packer packer;

    // set the buffer 
    packer.set_buffer(size, data);
    
    // pack the number of packed cells
    packer << num_packed_cells;

    // pack up the layout size
    packer << layout_size;

    // pack up the layout
    for (const int *i = packed_layout->begin(); i != packed_layout->end(); i++)
	packer << *i;

    // pack up the angle
    packer << theta_degrees;

    // pack up the extents
    for (int i = 0; i < extent_size; i++)
	packer << extent_data[i];

    Ensure (packer.get_ptr() == data + size);

    // clean up some memory
    delete [] extent_data;

    // make a packed mesh
    SP_Pack packed_mesh(new RZWedge_Mesh::Pack(size, data));

    Ensure (packed_mesh->get_num_packed_cells() == num_packed_cells);
    return packed_mesh;
}

//---------------------------------------------------------------------------//
/*!
  
 * \brief Pack up the cell extents.

 */
void RZWedge_Mesh::pack_extents(const sf_int &current_new,
				char *data,
				int   size,
				int   num_packed) const
{
    Require (current_new.size() == layout.num_cells());
    Require (cell_xz_extents.size() == layout.num_cells());
    Require (data != 0);
    Require (size == num_packed * 4 * sizeof(double));
    
    // loop through cells and create the new cell extents
    vf_double extents(num_packed, sf_double(4));
    for (int ncell, cell = 0; cell < cell_xz_extents.size(); cell++)
    {
	// find the new cell
	ncell = current_new[cell];
	
	if (ncell > 0)
	{
	    Check (ncell <= num_packed);
	    Check (extents[ncell-1].size() == cell_xz_extents[cell].size());

	    // add the extents to the new cell extents
	    for (int i = 0; i < extents[ncell-1].size(); i++)
		extents[ncell-1][i] = cell_xz_extents[cell][i];
	}
    }

    // now pack up the extents
    rtt_dsxx::Packer packer;
    packer.set_buffer(size, data);
    for (int i = 0; i < extents.size(); i++)
	for (int j = 0; j < extents[i].size(); j++)
	    packer << extents[i][j];

    Ensure (packer.get_ptr() == size + data);
} 

//===========================================================================//
// RZWEDGE_MESH::PACK DEFINITIONS
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct a RZWedge_Mesh::Pack instance.  Once allocated mesh data is
 * given to the RZWedge_Mesh::Pack constructor in the form of a char*, the
 * Pack object owns it.  When the Pack object goes out of scope it will clean
 * up the memory.  In general, Pack objects are only created by calling the
 * RZWedge_Mesh::pack() function.

 * \param s size of char data stream
 * \param d pointer to char data stream

 */
RZWedge_Mesh::Pack::Pack(int s, char *d)
    : data(d),
      size(s)
{
    Require (size >= (3 * sizeof(int) + sizeof(double)));
    Require (data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.

 * Do copy construction while preserving memory.  This is not a reference
 * counted class so data is copied from one class to the other during
 * function calls and the like (wherever a copy constructor is called).

 */
RZWedge_Mesh::Pack::Pack(const Pack &rhs)
    : data(new char[rhs.size]),
      size(rhs.size)
{
    // fill up new data array
    for (int i = 0; i < size; i++)
	data[i] = rhs.data[i];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.

 * Cleans up memory when the Pack object goes out of scope.  Once allocated
 * pointers are given to the Pack object the Pack object takes control of
 * them.

 */
RZWedge_Mesh::Pack::~Pack()
{
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get number of cells in the packed mesh.
 */
int RZWedge_Mesh::Pack::get_num_packed_cells() const
{
    Require (size >= (3 * sizeof(int) + sizeof(double)));

    int   num_cells = 0;

    rtt_dsxx::Unpacker unpacker;
    unpacker.set_buffer(size, data);

    unpacker >> num_cells;
    Check (unpacker.get_ptr() == data + sizeof(int));
    
    Ensure (num_cells >= 0);
    return num_cells;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the RZWedge_Mesh.

 * Unpacks and returns a smart pointer to the new RZWedge_Mesh.

 * \return smart pointer to the unpacked mesh

 */
RZWedge_Mesh::SP_Mesh RZWedge_Mesh::Pack::unpack() const
{
    using rtt_dsxx::SP;

    Require (size >= (3 * sizeof(int) + sizeof(double)));

    // build an XYZ coordinate system
    SP<Coord_sys> coord(new XYZCoord_sys());

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    unpacker.set_buffer(size, data);

    // determine the number of packed cells
    int num_packed_cells = 0;
    unpacker >> num_packed_cells;
    Check (num_packed_cells >= 0);

    // UNPACK THE LAYOUT
    int layout_size  = 0;
    unpacker >> layout_size;

    // don't need to reclaim this memory because we are giving it to the
    // layout packer
    int *layout_data = new int[layout_size];
    for (int *i = layout_data; i != layout_data + layout_size; i++)
	unpacker >> *i;

    AMR_Layout::Pack packed_layout(layout_size, layout_data);
    SP<AMR_Layout> layout = packed_layout.unpack();
    Check (layout->num_cells() == num_packed_cells);

    // GET THETA (degrees)
    double theta = 0;
    unpacker >> theta;
    Check (theta > 0);

    // UNPACK THE EXTENTS
    int       num_extents = num_packed_cells * 4;
    sf_double extent_data(num_extents, 0.0);
    for (int i = 0; i < num_extents; i++)
	unpacker >> extent_data[i];

    // build the new cell extents
    int       cectr = 0;
    vf_double cell_extents(num_packed_cells, sf_double(4));
    for (int i = 0; i < cell_extents.size(); i++)
	for (int j = 0; j < cell_extents[i].size(); j++)
	    cell_extents[i][j] = extent_data[cectr++];
    Check (cectr++ == num_extents);

    Ensure (unpacker.get_ptr() == size + data); 
    
    // build the new mesh
    SP<RZWedge_Mesh> mesh(new RZWedge_Mesh(coord, *layout, cell_extents,
					   theta, true));

    Ensure (mesh->num_cells() == num_packed_cells);
    Ensure (mesh->get_spatial_dimension() == coord->get_dim());
    Ensure (!mesh->full_Mesh());
    Ensure (mesh->get_total_volume() > 0.0);

    return mesh;
}

} // end namespace rtt_mc

#endif                          // rtt_mc_RZWedge_Mesh_cc

//---------------------------------------------------------------------------//
//                        end of mc/RZWedge_Mesh.cc
//---------------------------------------------------------------------------//
