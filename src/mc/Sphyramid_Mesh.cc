//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphyramid_Mesh.cc
 * \author Jeffery Densmore (Stolen from RZWedge_Mesh.cc)
 * \date   Mon Nov 10 2:16:00 2003
 * \brief  Sphyramid_Mesh implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Sphyramid_Mesh_cc
#define rtt_mc_Sphyramid_Mesh_cc

#include "Sphyramid_Mesh.hh"
#include "XYZCoord_sys.hh"
//#include "Constants.hh"
#include "viz/Ensight_Translator.hh"
#include "ds++/Packing_Utils.hh"
//#include <iomanip>
#include <typeinfo>
#include <algorithm>
#include <cmath>

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*! 
 * \brief Sphyramid_Mesh constructor.

 * The constructor requires complete arguments of the constituent data
 * necessary to build a Sphyramid_Mesh.  It is expected that an appropriate
 * builder class will message general data through an interface to build the 
 * specific data structures needed by Sphyramid_Mesh

 * \param coord_ coordinate system smart pointer

 * \param layout_ Layout giving cell connectivity information

 * \param cell_x_extents_ the x coordinates of each cell in the following
 * form: [cell][low x]; [cell][hi x];

 * \param beta_radians_ "angle" of Sphyramid in radians (not angle of
 * spherical
 
 * \param submesh_ [=false] boolean indicator; true if this is a submesh
 * cone)

*/

Sphyramid_Mesh::Sphyramid_Mesh(SP_Coord coord_, Layout &layout_, 
			       vf_double & cell_x_extents_, 
			       double beta_radians_, bool submesh_)
    :coord(coord_),
     layout(layout_),
     cell_x_extents(cell_x_extents_),
     beta_radians(beta_radians_),
     submesh(submesh_)
    

{
    using global::pi;

    // Check coordinate system class
    Require (this->coord);
    Require (this->coord->get_dim() == 3);
    Require (this->coord->get_Coord() == std_string("xyz"));
    
    // Make sure that the Sphyramid angle is positive and not too big
    Require (this->beta_radians > 0.0);
    Require (this->beta_radians <= pi/4.);

    // check that the cell-extents vector has num_cells elements
    // and weakly check that each cell-element has 2 extent-elements
    Check (this->cell_x_extents.size()    == this->layout.num_cells());
    Check (this->cell_x_extents[0].size() == 2);
    Check (this->cell_x_extents[this->layout.num_cells()-1].size() == 2);

    // precalculate heavily used angle quatities
    calc_angle_data();

    // calculate and assign the on-processor total volume
    calc_total_volume();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATIONS
//---------------------------------------------------------------------------//
/*! 
 * \brief pack up cell extents
 * 
 * \param current_to_new mapping of cells to be packed
 * \param data pointer to data stream
 * \param size size in bytes of cell extents to be packed
 * \param num_packed number of cells to be packed
 *
 */
void Sphyramid_Mesh::pack_extents(const sf_int &current_to_new,
				  char *data, int size, int num_packed) const
{
    Require (current_to_new.size()       == this->layout.num_cells());
    Require (this->cell_x_extents.size() == this->layout.num_cells());
    Require (data != 0);
    Require (size == 2*num_packed*sizeof(double));

    // loop through cells and create the new cell extents
    vf_double extents(num_packed, sf_double(2));
    for (int ncell, cell = 0; cell < this->cell_x_extents.size(); cell++)
    {
	// find the new cell
	ncell = current_to_new[cell];

	if (ncell > 0)
	{
	    Check (ncell <= num_packed);
	    Check (extents[ncell-1].size() == 
		   this->cell_x_extents[cell].size());

	    // add the extents to the new cell extents
	    for (int i = 0; i < extents[ncell-1].size(); i++)
	    {
		extents[ncell-1][i] = this->cell_x_extents[cell][i];
	    }
	}
    }

    // now pack up the extents
    rtt_dsxx::Packer packer;
    packer.set_buffer(size, data);
    for (int i = 0; i < extents.size(); i++)
    {
	for (int j = 0; j < extents[i].size(); j++)
	{
	    packer << extents[i][j];
	}
    }

    Ensure (packer.get_ptr() == data+size);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Calculate wedge angle data once for use throughout calculation
 *
 * \param beta_radians "angle" of Sphyramid in degrees (not angle of spherical
 * cone)
 */
void Sphyramid_Mesh::calc_angle_data()
{
    using global::pi;
    using std::tan;
    using std::sin;
    using std::cos;

    Require (this->beta_radians > 0.0);
    Require (this->beta_radians <= pi/4.);
   
    // precalculate heavily used trig functions
    this->tan_beta = std::tan(this->beta_radians);
    this->sin_beta = std::sin(this->beta_radians);
    this->cos_beta = std::cos(this->beta_radians);

    Ensure (this->tan_beta > 0.0);
    Ensure (this->sin_beta > 0.0);
    Ensure (this->cos_beta > 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and set the total (on-processor) volume of the mesh.
 *
 * \return
 */
void Sphyramid_Mesh::calc_total_volume()
{
    Require (num_cells()>0);

    // initialize private data member
    this->total_volume = 0.0;
    
    // sum local cell volumes
    for(int cell = 1; cell <= num_cells(); cell++)
	this->total_volume += volume(cell);

    Ensure (this->total_volume > 0.0);
}
//---------------------------------------------------------------------------//
// PUBLIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a position vector is within a cell.

 * \param cell cell index

 * \param r sf_double of position (3 dim vector)

 * \return true if r is in cell; false otherwise
 */
bool Sphyramid_Mesh::in_cell(int cell, const sf_double &r) const
{
    using rtt_mc::global::soft_equiv;

    Require (r.size() == 3);
    Require (cell > 0);
    Require (cell <= this->layout.num_cells());
    
    // first check x dimension
    if ((r[0] < this->cell_x_extents[cell-1][0] &&
	 !soft_equiv(r[0], this->cell_x_extents[cell-1][0]))
	||
	(r[0] > this->cell_x_extents[cell-1][1] &&
	 !soft_equiv(r[0], this->cell_x_extents[cell-1][1])))
	return false;

    // check y dimension
    if ((r[1] < -(r[0]*this->tan_beta) &&
	 !soft_equiv(r[1], -(r[0]*this->tan_beta)))
	||
	(r[1] > (r[0]*this->tan_beta) &&
	 !soft_equiv(r[1], -(r[0]*this->tan_beta))))
	return false;

    // check z dimension (cell is symmetric)
    if ((r[2] < -(r[0]*this->tan_beta) &&
	 !soft_equiv(r[2], -(r[0]*this->tan_beta)))
	||
	(r[2] > (r[0]*this->tan_beta) &&
	 !soft_equiv(r[2], -(r[0]*this->tan_beta))))
	return false;

    return true;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Return cell type for each cell in the Sphyramid_Mesh
 *
 * \return vector of strings cotaining cell type for each cell
 */
Sphyramid_Mesh::sf_int Sphyramid_Mesh::get_cell_types() const
{
    using std::fill;

    sf_int cell_type(this->layout.num_cells());

    // all cells in a Sphyramid_Mesh are general, 8-node hexahedrons
    fill(cell_type.begin(),cell_type.end(),rtt_viz::eight_node_hexahedron);

    return cell_type;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the coordinates of all nodes in the mesh.
 * 
 * For each cell in the Sphyramid_Mesh, the coordinates of all eight nodes are
 * returned. Thus, all interior nodes are replicated twice.  This approach,
 * although memory-inefficient, is easier to code especially considering that
 * the Sphyramid_Mesh does not already have the raw vertex data.
 *
 * \return vector of node coordinates for entire mesh
 *
 */
Sphyramid_Mesh::vf_double Sphyramid_Mesh::get_point_coord() const
{

    // number of vertices is always 8; always 3D
   const int num_verts_cell = 8;
   int vert_index;
   Check (this->coord->get_dim() == 3);

    // weakly check the validity of num_cells()
   Check(num_cells() >0 );

    //initialize the return vector
   vf_double return_coord(num_cells()*num_verts_cell);

    // for each cell, get vertices and reverse the columns and rows
   for (int cell = 1; cell <= num_cells(); cell++)
   {
	// get the cell's verices[dim][8]
        vf_double cell_verts = get_vertices(cell);

	// check validity of cell vertices vector
       Check (cell_verts.size() == this->coord->get_dim());

	// loop over all 8 nodes for this cell
       for (int node = 0; node < num_verts_cell; node++)
       {
	   //calculate running vertex index
	   vert_index = (cell-1)*num_verts_cell+node;

	   //resize each vertice's return_coord to num dimensions
	   return_coord[vert_index].resize(coord->get_dim());

	   // re-assign point coordinates to return vector
	   for (int dim = 0; dim < this->coord->get_dim(); dim++)
	   {
	       Check (cell_verts[dim].size() == num_verts_cell);
	       return_coord[vert_index][dim] = cell_verts[dim][node];
	    }
	}
    }

    // return the coordinate for all nodes of all cells
   return return_coord;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Find cell in Sphyramid_Mesh given position
 *
 * \param r position
 * 
 * \return cell containing position r, else -1 if r not in any cell
 */
int Sphyramid_Mesh::get_cell(const sf_double &r) const
{
    // This is some algorithm that Todd and Tom made up.
    // I'm including it so I don't get yelled at.

    using std::fabs;
    Require (!this->submesh);

    // initialize cell, flag indicating that cell was located
    int located_cell = 0;
    bool found       = false;

    // find cell that contains the x position
    while (!found)
    {
	located_cell += 1;
	double lox = get_low_x(located_cell);
	double hix = get_high_x(located_cell);
	
	if (r[0] >= lox && r[0] <= hix)
	    found=true;

	if (located_cell == num_cells() && !found)
	{
	    located_cell = -1;
	    found        = true;
	}
    }

    // check that the y- and z-position is with the cell
    if(located_cell > 0)
    {
        Check (fabs(r[1]) <= r[0]*this->tan_beta);
	Check (fabs(r[2]) <= r[0]*this->tan_beta);
    }

    return located_cell;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief calculate minimum distance to cell boundary
 * 
 * \param r position
 * \param omega direction
 * \param cell cell number
 * \param face index of boundary (return value)
 * \return mininum distance to boundary
 */
double Sphyramid_Mesh::get_db(const sf_double &r, const sf_double &omega,
			      int cell, int &face) const
{
    using global::dot;
    using global::soft_equiv;
    using global::huge;

    Require (r.size()     == 3);
    Require (omega.size() == 3);
    Require (soft_equiv(dot(omega,omega), 1.0, 1.0E-4));
    Require (in_cell(cell, r));

    // set up 6 dists-to-bndry, initialize to huge value
    sf_double distance(6,huge);

    // low x face (x-x_low=0)
    if (omega[0] < 0.0)
    {
	distance[0] = (get_low_x(cell)-r[0])/omega[0];
    }

    // high x face (x-x-high=0)
    if (omega[0] > 0.0)
    {
	distance[1] = (get_high_x(cell)-r[0])/omega[0];
    }

    // low y face (y+x*tan_beta=0)
    if (dot(omega,get_normal(cell,3)) > 0.0)
    {
	distance[2] = -(r[1]+r[0]*this->tan_beta)/
	    (omega[1]+omega[0]*this->tan_beta);
    }
    
    // high y face (y-x*tan_beta=0)
    if (dot(omega,get_normal(cell,4)) > 0.0)
    {
	distance[3] = (r[1]-r[0]*this->tan_beta)/
	    (omega[0]*this->tan_beta-omega[1]);
    }

    // low z face (z+x*tan_beta=0)
    if (dot(omega,get_normal(cell,5)) > 0.0)
    {
	distance[4] = -(r[2]+r[0]*this->tan_beta)/
	    (omega[2]+omega[0]*this->tan_beta);
    }
    
    // high z face (z-x*tan_beta=0)
    if (dot(omega,get_normal(cell,6)) > 0.0)
    {
	distance[5] = (r[2]-r[0]*this->tan_beta)/
	    (omega[0]*this->tan_beta-omega[2]);
    }


    // find face index and value of minimum distance
    double min_dist = huge;
    int face_index  = 0;

    for (int f = 0; f < 6; f++)
    {
	if (distance[f] < min_dist)
	{
	    min_dist  = distance[f];
	    face_index = f;
	}
    }
    
    // set face index
    Ensure (face_index >= 0);
    Ensure (face_index <  6);

    face = face_index+1;

    //return output
    Ensure (min_dist >= 0.0);
    Ensure (face > 0);
    Ensure (face <=6);

    return min_dist;
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the vertices for a given Sphyramid_Mesh cell.
 *
 * During normal use of the Sphyramid_Mesh, the explicit cell vertices are not
 * required. However, they are needed for graphics dumps.  The Sphyramid_Mesh
 * is always in 3-D XYZ geometry and always has eight vertices (the cell at
 * the radial center has four coincident vertices).
 * 
 * \param cell Sphyramid_Mesh (global) cell number
 *
 * \return vertices - the coordinates of the cell's eight vertices
 */
Sphyramid_Mesh::vf_double Sphyramid_Mesh::get_vertices(int cell) const
{
    Require (cell > 0);
    Require (cell <= num_cells());
    Check   (this->coord->get_dim() == 3);

    const int num_verts_face = 4;
    const int num_verts_cell = 8;
    const int loz_face       = 5;
    const int hiz_face       = 6;

    // get the vertices for the low z face of the cell
    vf_double cell_vertices = get_vertices(cell, loz_face);

    // get the vertices for the high z face of the cell
    vf_double hiz_face_vertices = get_vertices(cell,hiz_face);

    Check (cell_vertices.size()    == this->coord->get_dim());
    Check (cell_vertices[0].size() == num_verts_face);
    Check (cell_vertices[1].size() == num_verts_face);
    Check (cell_vertices[2].size() == num_verts_face);

    Check (hiz_face_vertices.size()    == this->coord->get_dim());
    Check (hiz_face_vertices[0].size() == num_verts_face);
    Check (hiz_face_vertices[1].size() == num_verts_face);
    Check (hiz_face_vertices[2].size() == num_verts_face);

    // add the hiz face vertices to the loz face vertices
    for (int v = 0; v < this->coord->get_dim(); v++)
    {
	cell_vertices[v].insert(cell_vertices[v].end(),
			hiz_face_vertices[v].begin(),
			hiz_face_vertices[v].end());
    }

    // make sure each dimension has 8 entries -- one per vertex
    Ensure (cell_vertices[0].size() == num_verts_cell);
    Ensure (cell_vertices[1].size() == num_verts_cell);
    Ensure (cell_vertices[2].size() == num_verts_cell);

    // return the cell vertices
    return cell_vertices;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the vertices for a given face of a Sphyramid_Mesh cell. 
 * 
 * During normal use of the Sphyramid_Mesh, the explicit cell vertices are not
 * required. However, they are needed for graphics dumps.  The Sphyramid_Mesh
 * is always in 3-D XYZ geometry and each face always has four vertices (the
 * cell at the radial center has four coincident vertices).
 * 
 * \param cell Sphyramid_Mesh (global) cell number
 * \param face face of cell
 *
 * \return vertices of the coordinates of the four vertices defining a face
 */
Sphyramid_Mesh::vf_double Sphyramid_Mesh::get_vertices(int cell, int face) const
{
    Require (face > 0);
    Require (face <= 6);
    Require (cell > 0);
    Require (cell <= this->layout.num_cells());
    Check (this->coord->get_dim() == 3);

    vf_double face_vertices(this->coord->get_dim());
    sf_double single_vert(this->coord->get_dim(),0.0);

    double lox     = get_low_x(cell);
    double hix     = get_high_x(cell);
    double small_y = lox*this->tan_beta;
    double large_y = hix*this->tan_beta;
    double small_z = lox*this->tan_beta;
    double large_z = hix*this->tan_beta;

    // low x face or high x face
    if (face == 1 || face == 2)
    {
	double y_variable;
	double z_variable;

	if (face == 1)
	{
	    single_vert[0]=lox;
	    y_variable=small_y;
	    z_variable=small_z;
	}
	else if (face == 2)
	{
	    single_vert[0] = hix;
	    y_variable = large_y;
	    z_variable= large_z;
	}
	
	single_vert[1] = -y_variable;
	single_vert[2] = -z_variable;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[1] = y_variable;
	for(int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[2]= z_variable;
	for(int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[1]=-y_variable;
	for(int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
    }

    // low y face or high y face
    else if (face == 3 || face == 4)
    {
	double y_plusminus;
	if (face == 3)
	{
	    y_plusminus = -1.0;
	}
	else if (face == 4)
	{
	    y_plusminus = 1.0;
	}
	single_vert[0] = lox;
	single_vert[1] = y_plusminus*small_y;
	single_vert[2] = -small_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[0] = hix;
	single_vert[1] = y_plusminus*large_y;
	single_vert[2] = -large_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[0] = hix;
	single_vert[1] = y_plusminus*large_y;
	single_vert[2] = large_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);	
	}
	single_vert[0] = lox;
	single_vert[1] = y_plusminus*small_y;
	single_vert[2] = small_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);	
	}
    }

    // low z face or high z face
    
    else if (face == 5 || face == 6)
    {
	double z_plusminus;
	if (face == 5)
	{
	    z_plusminus = -1.0;
	}
	else if (face == 6)
	{
	    z_plusminus = 1.0;
	}
	single_vert[0] = lox;
	single_vert[1] = -small_y;
	single_vert[2] = z_plusminus*small_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[0] = hix;
	single_vert[1] = -large_y;
	single_vert[2] = z_plusminus*large_z;
	for (int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);
	}
	single_vert[0] = hix;
	single_vert[1] = large_y;
	single_vert[2] =z_plusminus*large_z;
	for(int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);	
	}
	single_vert[0] = lox;
	single_vert[1] = small_y;
	single_vert[2] = z_plusminus*small_z;
	for(int i = 0; i < this->coord->get_dim(); i++)
	{
	    face_vertices[i].push_back(single_vert[i]);	
	}
    }

    // check that there are 3 dimensions a 4 vertices
    Ensure (face_vertices.size()    == this->coord->get_dim());
    Ensure (face_vertices[0].size() == 4);
    Ensure (face_vertices[1].size() == 4);
    Ensure (face_vertices[2].size() == 4);
 
    // return the four vertices for the face
    return face_vertices;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief return the face number for a given cell boundary (independent of cell)
 * 
 * \param boundary string boundary description
 * \param cell cell number
 *
 * \return face number
 */

int Sphyramid_Mesh::get_bndface(std_string boundary, int cell) const
{
    //return value
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
	face =6;
    else
	Insist (0,"Unknown boundary string used!");

    // return the face number given the cell boundary
    return face;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief return a list of cells along a specified boundary
 * 
 * \param boundary string boundary name
 * \return vector of cell numbers
 */

Sphyramid_Mesh::sf_int Sphyramid_Mesh::get_surcells(std::string boundary) const
{

    Require(!this->submesh);
    Check (this->coord->get_dim() == 3);

    // make return vector containing a list of cells along specified boundary
    sf_int return_list;

    // verify assumption that cell 1 is the low x cell
    Insist (layout(1,1) == 1, "Cell 1 is not reflective on low x face!");
    Insist (layout(1,3) == 1, "Cell 1 is not reflective on low y face!");
    Insist (layout(1,4) == 1, "Cell 1 is not reflective on high y face!");
    Insist (layout(1,5) == 1, "Cell 1 is not reflective on low z face!");
    Insist (layout(1,6) == 1, "Cell 1 is not reflective on high z face!");
    
    // calculate the cells along the...
    // ... the high r boundary
    if ((boundary == "hix") || (boundary == "hir"))
    {
	// first, work our way along mesh to the high x cell
	int start_cell = 1;
	while (layout(start_cell,2) != start_cell && layout(start_cell,2) !=0)
	{
	    //get the next cell
	    int next_cell=layout(start_cell,2);

	    // ridiculous, unnecessary checks on the next cell
	    Check (next_cell != 0);
	    Check (next_cell != start_cell);
	    Check (next_cell  > 0);
	    Check (next_cell <= this->layout.num_cells());

	    //check that the cell has reflecting sides
	    Check (this->layout(next_cell,3) == next_cell);
	    Check (this->layout(next_cell,4) == next_cell);
	    Check (this->layout(next_cell,5) == next_cell);
	    Check (this->layout(next_cell,6) == next_cell);

	    // the next cell is okay
	    start_cell = next_cell;
	}
	
	// we have the high x cell;
	int current_cell=start_cell;
	return_list.push_back(current_cell);

	// check the size of the surface cell list
	Ensure (return_list.size() == 1);
    }
    else{
	Insist(0,
	       "Unkown or invalid (lor/lox,loy/hiy,loz/hiz) "
	       "surf in get_surcells!");
    }

    // return vector
    return return_list;
}
//---------------------------------------------------------------------------//
/*! * \brief check that a user/host-defined set of surface source cells
 *  actually resides on the surface of the system (requires a vacuum bnd).
 * 
 * \param ss_face boundary face
 * \param ss_list reference to a vector of cells
 * \return true if cells comprise a valide surface source
 */
bool Sphyramid_Mesh::check_defined_surcells(const std_string ss_face,
					 const sf_int &ss_list) const
{
    // a weak check on the number of surface cells
    Check (ss_list.size() <= num_cells());

    for (int ss_indx = 0; ss_indx < ss_list.size(); ss_indx++)
    {
	// convert face on which ss resides from string to int.
	// despite its args, get_bndface actually has no cell dependence
	int ss_face_num = get_bndface(ss_face, ss_list[ss_indx]);

	// get_bnd condition on ss face; had better be vacuum (0)
	int bc = this->layout(ss_list[ss_indx], ss_face_num);
	if (bc != 0)
	    return false;
    }

    return true;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Sample a position on the surface of a sphere inside a cell
 * 
 * This function is used by the Random Walk procedure to move a particle to
 * the edge of the RW sphere.  The sphere is always inside, and never on, the
 * cell boundaries.
 *
 * \param cell cell index
 * \param origin sphere origin
 * \param radius sphere radius
 * \param random random number object
 * \return a pair of vector<double> where the first element is the position
 * on the surface of the sphere, and the second element is the normal of the 
 * sphere at that direction.
 */
Sphyramid_Mesh::pair_sf_double Sphyramid_Mesh::sample_random_walk_sphere(
    int cell, const sf_double &origin, double radius, rng_Sprng &random) const
{
    using rtt_mc::global::dot;
    using std::make_pair;
    using rtt_mc::global::soft_equiv;
    
    Require (cell > 0);
    Require (cell <= this->layout.num_cells());
    Require (origin.size() == 3);
    Require (in_cell(cell, origin));

    // checks to make sure entire sphere is in cell
    Require (origin[0]+radius < get_high_x(cell));
    Require (origin[0]-radius > get_low_x(cell));

    // get initial position and direction of track
    sf_double     r = origin;
    sf_double omega = this->coord->sample_isotropic_dir(random);

    // distance we have to track
    double track = radius;

    // track until we have gone the radial distance
    double d_bnd;
    int face;
    double factor;
    sf_double normal;
    while (track > 0.0)
    {
	// determine shortest distance to boundary
	d_bnd = get_db(r, omega, cell, face);

	// process a reflection (y or z face)
	if (d_bnd < track)
	{

	    // make sure random walk sphere is truncated only by y or z faces
	    Check (face == 3 || face == 4 || face == 5 || face == 6);

	    // stream to the face
	    r[0] = r[0]+d_bnd*omega[0];
	    r[1] = r[1]+d_bnd*omega[1];
	    r[2] = r[2]+d_bnd*omega[2];

	    // adjust the remaining track_length
	    track -= d_bnd;
	    Check (track >= 0.0);

	    // calculate face normal
	    normal = get_normal(cell, face);
	    Check (normal.size() == 3);

	    // specularly reflect angle
	    factor    = dot(omega, normal);
	    omega[0] -= 2.0*factor*normal[0];
	    omega[1] -= 2.0*factor*normal[1];
	    omega[2] -= 2.0*factor*normal[2];
	}
	// stream until finished
	else
	{
	    r[0] = r[0]+track*omega[0];
	    r[1] = r[1]+track*omega[1];
	    r[2] = r[2]+track*omega[2];

	    // we are finished tracking
	    track = 0.0;
	}
    }
    // assign return position and normal
    pair_sf_double pos_and_norm = make_pair(r, omega);

    Ensure (in_cell(cell, pos_and_norm.first));
    Ensure (soft_equiv(dot(pos_and_norm.second, pos_and_norm.second), 1.0));

    return pos_and_norm;

}
//---------------------------------------------------------------------------//
// Mesh Packing Interface
//---------------------------------------------------------------------------//
/*! 
 * \brief Pack up a mesh into a Pack struct for communication and
 * persistence.
 *
 * The cell list provides the cells to pack up.  It also is a map from a full
 * mesh to a spatially decomposed mesh.  Thus, the packed mesh will only
 * contain cells in the cell list with the provided mappings.
 *
 * To produce an exact copy, call pack without any arguments.
 *
 * \param current_to_new list of cells to include in the packed mesh, set
 * this to NULL to produce an exact copy
 *
 * \return packed mesh
 */
Sphyramid_Mesh::SP_Pack Sphyramid_Mesh::pack(const sf_int &current_to_new)
    const
{
    Require (current_to_new.size() == this->layout.num_cells() ||
	     current_to_new.size() == 0);

    // determine wheter this is exact replication or sub-packing
    sf_int current_to_new_replicate;
    bool replicate;
    
    if (current_to_new.size() == 0)
    {
	replicate = true;
	current_to_new_replicate.resize(this->layout.num_cells());
	
	for (int cell = 1; cell <= this->layout.num_cells(); cell++)
	{
	    current_to_new_replicate[cell-1] = cell;
	}
    }
    else
    {
	replicate = false;
    }
    
    // the coordinate system is always XYZ
    Check (typeid(*(this->coord)) == typeid(XYZCoord_sys));

    // packup layout data
    Layout::SP_Pack packed_layout;
    int layout_size;
    int num_packed_cells;

    // packed extents data
    int   extent_size = 0;
    char *extent_data = 0;

    // packup the layout
   
    if (replicate)
    {
	// packup the layout
	packed_layout = this->layout.pack(current_to_new_replicate);
	layout_size = packed_layout->get_size();
	Check (layout_size >= 1);
 
	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells == this->layout.num_cells());

	// calculate extent size
	extent_size = 2*num_packed_cells*sizeof(double);
	extent_data = new char[extent_size];
	

	// pack extents
	pack_extents(current_to_new_replicate, extent_data, extent_size,
		     num_packed_cells);

    }
    else
    {
	// packup the layout
	packed_layout = this->layout.pack(current_to_new);
	layout_size   = packed_layout->get_size();
	Check (layout_size >= 1);

	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells <= this->layout.num_cells());

	// calculate extent size
	extent_size = 2*num_packed_cells*sizeof(double);
	extent_data = new char[extent_size];

	// pack extents
	pack_extents(current_to_new, extent_data, extent_size, 
		     num_packed_cells);
    }
    Check (num_packed_cells >= 0);
    Check (num_packed_cells <= this->layout.num_cells());
    Check (extent_size == 2*num_packed_cells*sizeof(double));
    Check (extent_data != 0);

    // now pack up the mesh

    // ints (1-size of packed_layout, 1-num_packed_cells, packed layout)
    int total_ints = (2+layout_size)*sizeof(int);

    // doubles (1-beta_radians)
    int total_doubles = 1*sizeof(double);

    // chars
    int total_chars = extent_size;

    // allocate space
    int size = total_ints + total_doubles + total_chars;
    char *data = new char[size];
    
    // pack up the mesh
    rtt_dsxx::Packer packer;

    // set the buffer
    packer.set_buffer(size, data);
   
    // pack the number of packed cells
    packer << num_packed_cells;

    // pack up the layout size
    packer << layout_size;

    // pack up the layout
    for (const int *i = packed_layout->begin(); i != packed_layout->end(); i++)
    {
	packer <<  *i;
    }

    // pack up the angle
    packer << beta_radians;

    // pack up the extents
    for (int i = 0; i < extent_size; i++)
    {
	packer << extent_data[i];
    }

    Ensure (packer.get_ptr() == data+size);

    // clean up some memory
    delete [] extent_data;

    // make a packed mesh
    SP_Pack packed_mesh(new Sphyramid_Mesh::Pack(size, data));
    
    Ensure (packed_mesh->get_num_packed_cells() == num_packed_cells);

    return packed_mesh;
}

//---------------------------------------------------------------------------//
// Overloaded == operator
//---------------------------------------------------------------------------//
bool Sphyramid_Mesh::operator==(const Sphyramid_Mesh &rhs) const
{
    using rtt_mc::global::soft_equiv;

    // check to see that the Layouts are equal
    if (this->layout != rhs.layout) return false;

    // check the x extents
    for (int cell = 1; cell <= this->layout.num_cells(); cell++)
    {
	if(!soft_equiv(get_low_x(cell),rhs.get_low_x(cell))) return false;

	if(!soft_equiv(get_high_x(cell),rhs.get_high_x(cell))) return false;
    }

    // check the beta angle
    if(!soft_equiv(this->beta_radians,rhs.beta_radians)) return false;

    // don't bother with coord system, all Sphyramid Meshs are XYZ

    // if we haven't returned, then the two meshes must be equal

    return true;
}
//===========================================================================//
// SPHYRAMID_MESH::PACK DEFINITIONS
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct a Sphyramid_Mesh::Pack instance.  Once allocated mesh data is
 * given to Sphyramid_Mesh::Pack constructor in the form of a char*, the
 * Pack object owns it.  When the Pack object goes out of scope it will clean
 * up the memory.  In general, Pack objects are only created by calling the
 * Sphyramid_Mesh::pack() function.

 * \param s size of char data stream
 * \param d pointer to char data stream

 */
Sphyramid_Mesh::Pack::Pack(int s, char *d)
    :data(d),
     size(s)
{
    Require (this->size >= (3*sizeof(int)+sizeof(double)));
    Require (this->data);
}
//---------------------------------------------------------------------------//
/*! 
 * \brief copy constructor
 * 
 * do copy consruction while preserving memory.  This is not a reference
 * counted class so data is copied from one object to the other during
 * function calls, etc.
 * 
 */
Sphyramid_Mesh::Pack::Pack(const Pack &rhs)
    :data(new char[rhs.size]),
     size(rhs.size)
{

    // fill up new data array
    for (int i = 0; i < size; i++)
    {
	this->data[i] = rhs.data[i];
    }

}
//---------------------------------------------------------------------------//
/*! 
 * \brief Destructor
 * 
 * Cleans up memory when the Pack object goes out of scope.  Once allocated
 * pointers are given to the Pack object the Pack object takes control of
 * them.
 */
Sphyramid_Mesh::Pack::~Pack()
{
    delete [] this->data;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief get number of cells in the packed mesh
 * 
 */
int Sphyramid_Mesh::Pack::get_num_packed_cells() const
{
    Require (this->size >= (3 *sizeof(int)+sizeof(double)));

    int num_cells = 0;

    rtt_dsxx::Unpacker unpacker;
    unpacker.set_buffer(this->size, this->data);
   
    unpacker >> num_cells;
      

    Check (unpacker.get_ptr() == data + sizeof(int));

    Ensure (num_cells >= 0);
    
    return num_cells;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Unpack the Sphyramid_Mesh
 * 
 * Unpacks and returns a smart pointer to the new Sphyramid_Mesh
 *
 * \return smart pointer to the unpacked mesh
 */
Sphyramid_Mesh::SP_Mesh Sphyramid_Mesh::Pack::unpack() const
{
    using rtt_dsxx::SP;

    Require (this->size >= (3*sizeof(int)+sizeof(double)));

    // build an XYZ coordinate system
    SP<Coord_sys> coord(new XYZCoord_sys());

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    unpacker.set_buffer(this->size, this->data);

    // determine the number of packed cells
    int num_packed_cells = 0;
    unpacker >> num_packed_cells;
    Check (num_packed_cells >= 0);   

    // unpack layout
    int layout_size = 0;
    unpacker >> layout_size;

    // don't need to reclaim this memory becuase we are giving it to the
    // layout packer
    int *layout_data = new int[layout_size];
    for (int *i = layout_data; i != layout_data+layout_size; i++)
    {
	unpacker >> *i;
    }

    Layout::Pack packed_layout(layout_size, layout_data);
    SP<Layout> layout = packed_layout.unpack();
    Check (layout->num_cells() == num_packed_cells);

    // get beta (radians)
    double beta_radians = 0;
    unpacker >> beta_radians;
    Check (beta_radians > 0);

    // unpack the extents
    int num_extents = num_packed_cells*2;
    sf_double extent_data(num_extents, 0.0);
    for (int i = 0; i < num_extents; i++)
    {
	unpacker >> extent_data[i];
    }

    // build the new cell extents
    int cectr = 0;
    vf_double cell_extents(num_packed_cells, sf_double(2));
    for (int i = 0; i < cell_extents.size(); i++)
    {
	for(int j = 0; j < cell_extents[i].size(); j++)
	{
	    cell_extents[i][j] = extent_data[cectr++];
	}
    }
    Check (cectr++ == num_extents);
    Check (unpacker.get_ptr() == this->size + this->data);

    // build the new mesh
    SP<Sphyramid_Mesh> mesh(new Sphyramid_Mesh(coord, *layout, cell_extents,
					       beta_radians,true));

    Ensure (mesh->num_cells()             == num_packed_cells);
    Ensure (mesh->get_spatial_dimension() == coord->get_dim());
    Ensure (!mesh->full_Mesh());
    Ensure (mesh->get_total_volume()      >  0.0);

    return mesh;
}

} // end namespace rtt_mc

#endif                        // rtt_mc_Sphyramid_Mesh_cc

//---------------------------------------------------------------------------//
//                 end of Sphyramid_Mesh.cc
//---------------------------------------------------------------------------//
