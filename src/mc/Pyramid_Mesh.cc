//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Pyramid_Mesh.cc
 * \author Jeffery Densmore
 * \date   Mon Oct  6 09:15:12 2003
 * \brief  Pyramid_Mesh implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Pyramid_Mesh_cc
#define rtt_mc_Pyramid_Mesch_cc

#include "Pyramid_Mesh.hh"
#include "Constants.hh"
#include "Math.hh"
#include "viz/Ensight_Translator.hh"

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*! 
 * \brief Pyramid_Mesh constructor.

 * The constructor requires complete arguments of the constituent data
 * necessary to build a Pyramid_Mesh.  It is expected that an appropriate
 * builder class will message general data through an interface to build the 
 * specific data structures needed by Pyramid_Mesh

 * \param coord_ coordinate system smart pointer

 * \param layout_ AMR_layout giving cell connectivity information

 * \param cell_x_extents_ the x coordinates of each cell in the following
 * form: [cell][low x]; [cell][hi x];

 * \param beta_degrees_ "angle" of pyramid in degrees (not angle of spherical
 * cone)

*/

Pyramid_Mesh::Pyramid_Mesh(SP_Coord coord_,
			   AMR_Layout &layout_,
			   vf_double & cell_x_extents_,
			   double beta_degrees_)
    :coord(coord_),
     layout(layout_),
     cell_x_extents(cell_x_extents_),
     beta_degrees(beta_degrees_)
{
    // Check coordinate system class
    Require (coord);
    Require (coord->get_dim() == 3);
    Require (coord->get_Coord() == std:string("xyz"));

    // Make sure that the pyramid angle is postive and not obtuse
    Require ((beta_degrees>0.0) && (beta_degrees <= 90.0));

    // check that the cell-extents vector has num_cells elements
    // and weakly check that each cell-element has 2 extent-elements
    Check (cell_x_extents.size() == layout.num.cells());
    Check (cell_x_extents[0].size() == 2);
    Check (cell_x_extents[layout.num_cells()-1].size() ==2);

    // precalculate heavily used angle quatities
    calc_angle_data(beta_degrees);

    // calculate and assign the on-processor total volume
    calc_total_volume()
}

 
//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate wedge angle data once for use throughout calculation
 *
 * \param beta_degrees_ "angle" of pyramid in degrees (not angle of spherical
 * cone)
 */
void Pyramid_Mesh::calc_angle_data(const double beta_degrees)
{
    Require ((beta_degrees >0.0) && (beta_degrees <= 90.0));
    
    beta_radians=beta_degrees*rtt_mc::global::pi/180.0;
    Check (beta_radians>0.0);
    Check (beta_radians <=rtt_mc::global::pi/2.0);

    // precalculate heavily used trif functions
    tan_beta=std::tan(beta_radians);
    sin_beta=std::sin(beta_radians);
    cos_beta=std::cos(beta_radians);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate and set the total (on-processor) volume of the mesh.
 *
 * \return total volume of on-processor Pyramid cells
 */
void Pyramid_Mesh::calc_total_volume()
{
    Require (num_cells()>0);

    // initialize private data member
    total_volume=0.0;
    
    // sum local cell volumes
    for(int cell=1; cell<=num_cells(); cell++)
	total_volume+=volume(cell);

    Ensure (total_volume >0.0);
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
bool Pyramid_Mesh::in_cell(int cell, const sf_double &r) const
{
    using rtt_mc::global::soft_equiv;

    Require (r.size() == 3);
    Require (cell >0 && cell <= layout.num_cells());
    
    // first check x dimension
    if ((r[0]<cell_x_extents[cell-1][0] &&
	 !soft_equiv(r[0],cell_x_extents[cell-1][0]))
	||
	(r[0]> cell_x_extents[cell-1][1] &&
	 !soft_equiv(r[0],cell_x_extents[cell-1][1])))
	return false;

    // check y dimension
    if ((r[1] < -(r[0]*tan_beta) &&
	 !soft_equiv(r[1], -(r[0]*tan_beta)))
	||
	(r[1] > (r[0]*tan_beta) &&
	 !soft_equiv(r[1], -(r[0]*tan_beta))))
	return false;

    // check z dimension (cell is symmetric)
    if ((r[2] < -(r[0]*tan_beta) &&
	 !soft_equiv(r[2], -(r[0]*tan_beta)))
	||
	(r[2] > (r[0]*tan_beta) &&
	 !soft_equiv(r[2], -(r[0]*tan_beta))))
	return false;

    return true;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR IMC
//---------------------------------------------------------------------------//
/*!
 * \Brief Calculate minimum distance to boundary in an Pyramid_Mesh cell
 *
 * \param r position
 * \param omega direction
 * \param cell cell containing position r
 *
 * \return get_db minimum distance to boundary
 * \return face (in arguments) face corresponding to min dist-to-bndry
 */
double Pyramid_Mesh::get_db(const sf_double &r, const sf_double *omega,
			    int cell, int &face) const
{
    using std::vector;
    using global::dot;
    Check (r.size() == 3);
    Check (omega.size() ==3);
    Check (global::soft_equiv(dot(omega,omega),1.0, 1.0e-5));

    // set up 6 dists-to-bndry, initialize to huge value
    // -- always 6 faces in an Pyramid_Mesh cell.
    vector<double> distance(6,global::Huge);

    // low x face (inner radial face) (internal index 0)
    if(omega[0]< 0.0)
	distance[0]=(get_low_x(cell)-r[0])/omega[0];

    // high x face (outer radial face) (internal index 1)
    if(omega[0]>0.0)
	distance[1]=(get_high_x(cell) - r[0]) / omega[0];
    
    // low y face (tan(beta)x+y=0) (internal index 2)
    if(dot(omega,get_normal(cell,3))>0.0)
	distance[2] = -(r[1]+r[0]*tan_beta)/
	    (omega[1]+omega[0]*tan_beta);

    // high y face (tan(beta)x-y=0) (internal index 3)
    if(dot(omega,get_normal(cell,4))>0.0)
	distance[3]=(r[1]-r[0]*tan_beta)/
	    (omega[0]*tan_beta-omega[1]);

    // low z face (tan(beta)x+z=0) (internal index 4)
    if(dot(omega,get_normal(cell,5))>0.0)
	distance[4] = -(r[2]+r[0]*tan_beta)/
	    (omega[2]+omega[0]*tan_beta);

    // high z face (tan(beta)x-z=0) (internal index 5)
    if(dot(omega,get_normal(cell,6))>0.0)
	distance[5]=(r[2]-r[0]*tan_beta)/
	    (omega[0]*tan_beta-omega[2]);

    // find face index and value of minimum(distance[face])
    double min_dist = global::huge;
    int face_index =0 ;
    for (int f=0; f<6; f++)
	if (distanc[f]<min_dist)
	{
	    min_dis = distance[f];
	    face_index=f;
	}

    // set face index (external: in [1,6])
    Ensure (face_index>=0 && face_index <6);
    face=face_index+1;

    // check that return quantities are within expected limits
    Ensure (min_dist >= 0.0);
    Ensure (face>0 && face<=6);

    // return minimum distance to boundary
    return min_dist;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Find cell in Pyramid_Mesh give position
 *
 * \param r position
 * 
 * \return cell cell containing position r
 */
int Pyramid_Mesh::get_cell(const sf_double &r) const
{
    // This is some algorithm that Todd and Tom made up.
    // I'm including it so I don't get yelled at.

    // initialize cell, flag indicating that cell was located
    in located_cell =0;
    bool found =false;

    // find cell that contains the x position
    while(!found)
    {
	located_cell+=1;
	double lox = get_low_x(located_cell);
	double hix = get_high_x(located_cell);
	
	if (r[0]>= lox && r[0] <=hix)
	    found=true;

	if (located_cell == num_cells() && !found)
	{
	    located_cell= -1;
	    found=true;
	}
    }

    // check that the y- and z-position is with the cell
    if(located_cell>0)
    {
	Check(std::fabs(r[1])<=r[0]*tan_beta);
	Check(std::fabs(r[2])<=r[0]*tan_beta);
    }

    return located_cell;
}
//---------------------------------------------------------------------------//
// return the face number for a given cell boundary (independent of cell)

int Pyramid_Mesh::get_bndface(std_string boundary, int cell) const
{
    //return value
    int face;

    if (boundary =="lox" || boundary == "lor")
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
// return a list of cells along a specified boundary

Pyramid_Mesh::sf_int Pyramid_Mesh::get_surcells(std::string boundary) const
{
    std::vector;

    Require(coord->get_dim() == 3);

    // make return vector containing a list of cells along specified boundary
    vector<int> return_list;

    // verify assumption that cell 1 is the low x cell
    Insist ((layout(1,1,1)==1) && (layout(1,3,1)==1) && (layaout(1,4,1)==1)
	    && (layout(1,5,1)==1) && (layout(1,6,1)==1),
	    "Cell 1 is not reflect on inside and sides!");
    
    // calculate the cellss along the...
    // ... the high r boundary
    if ((boundary == "hix") || (boundary == "hir"))
    {
	// first, work our way along mesh to the high x cell
	int start_cell = 1;
	while(layout(start_cell,2,1) != start_cell &&
	      layout(start_cell,2,1) !=0)
	{
	    //get the next cell
	    int next_cell=layout(start_cell,2,1);

	    // ridiculous, unnecessary checks on the next cell
	    Check (next_cell !=0 && next_cell != start_cell);
	    Check (next_cell >0 && next_cell <= layout.num_cells());

	    //check that the cell has reflecting sides
	    Check ((layout(next_cell,3,1) == next_cell) &&
		   (layout(next_cell,4,1) == next_cell) &&
		   (layout(next_cell,5,1) == next_cell) &&
		   (layout(next_cell,6,1) == next_cell));

	    // the next cell is okay
	    start_cell = next_cell;
	}
	
	// we have the high x cell;
	int current_cell=start_cell;
	return_list.push_back(current_cell);

	// check the size of the surface cell list
	Check (return_list.size() == 1);
    }
    else
	Insist(0,
	"Unkown or invalid (lor/lox,loy/hiy,loz/hiz) surf in get_surcells!");

    // return vector
    return return_list;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Sample a position on the surface of a sphere inside a cell.
 *
 * This function is required by the IMC_MT concept. This function is used
 * exclusively by the random walk procedure, wherein the sphere should always
 * be inside, and never on, the cell boundaries.
 *
 * As opposed to the OS_Mesh Implementation of this function, this function
 * actually tracks a distance equal to radius on the Pyramid_Mesh.  The
 * returned normal is actally the direction of the ray when it reaches that
 * distance
 * 
 * \param cell cell index
 * \param origin sphere origin
 * \param radius sphere radius
 * \param random random number object
 * 
 * \return a [pair of vector<double> where the first element of the pair the
 * position on the surface of the sphere and the second element of the pair
 * is the normal of the sphere at that position
 */
Pyramid_Mesh::pair_sf_double Pyramid_Mesh::sample_random_walk_sphere(
    int cell,
    const sf_double & &origin,
    double radius,
    rng_Sprng &random) const
{
    Require (cell >0);
    Require (cell <= layout.num_cells());
    Require (origin.size() == 3);
    Require (in_cell(cell,origin));

    // checks to make sure sphere is in cell in x dimension
    Require(origin[0]-radius>get_low_x(cell));
    Require(origin[0]+radius<get_high_x(cell));

    // get initial position and direction to track
    sf_double r = origin;
    sf_double omega=coord->sample_isotropic_dir(random);

    // distance we have to track
    double track = radius;

    // track_until we have gone the radial distance
    double d_bnd =0.0;
    int face = 0;
    double factor =0.0;
    sf_double normal;
    while (track >0.0)
    {
	// determine shortest distance to boundary
	d_bnd=get_db(r,omega,cell,face);

	// process a reflection on y or z face
	if(d_bnd< track)
	{

	    // make sure random walk sphere is not truncated by x face
	    Check (face == 3 || face==4 || face==5 || face==6);

	    // stream to face
	    r[0]=r[0]+d_bnd*omega[0];
	    r[1]=r[1]+d_bnd*omega[1];
	    r[2]=r[2]+d_bnd*omega[2];

	    // adjust the remaining track length
	    track-=d_bnd;
	    Check (track >= 0.0);

	    // calculate face normal
	    normal=get_normal(cell,face);
	    Check (normal.size() == 3);

	    // specularly reflect angle
	    factor=rtt_mc::global::dot(omega,normal);
	    omega[0]-=2.0*factor*normal[0];
	    omega[1]-=2.0*factor*normal[1];
	    omega[2]-=2.0*factor*normal[2];
	}

	// stream until we are finished

	else
	{
	    r[0]=r[0]+track*omega[0];
	    r[1]=r[1]+track*omega[1];
	    r[2]=r[2]+track*omega[2];

	    // we are finished tracking
	    track =0.0;
	}
    }

    // assign and return position and normal, the normal is the last
    // direction the ray had before reaching the specified track distance
    pair_sf_double pos_and_norm= std::make_pair(r,omega);

    // checks; we do not check that omega is properly normalize here because
    // it is checked in the preceding and ensuing transport steps
    Ensure (in_cell(cell, pos_and_norm.first));

    // return
    return pos_and_norm;
}
//---------------------------------------------------------------------------//
// check that a user/host-defined set of surface source cells actually
// resides on the surface of the system (requires a vacuum bnd).

bool Pyramid_Mesh::check_define_surcells(const std_string ss_face,
					 const sf_int &ss_list) const
{
    // a weak check on the number of surface cells
    Check (ss_list.size() <= num_cells());

    for (int ss_indx = 0; ss_indx<ss_list.size(); ss_indx++)
    {
	// convert face on which ss resides from string to in.
	// despite its args, get_bndface actually has no cell dependence
	int ss_face_num = get_bndface(ss_face,ss_list[ss_indx]);

	// get_bnd condition on ss face; had better be vacuum (0)
	int bc = layout(ss_list[ss_indx],ss_face_num,1);
	if (bc!=0)
	    return false;
    }

    return true;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the vertices for a given Pyramid_Mesh cell.
 *
 * During normal use of the Pyramid_Mesh, the explicit cell vertices are not
 * required. However, they are needed for graphics dumps.  The Pyramid_Mesh
 * is always in 3-D XYZ geometry and always has eight vertices (the cell at
 * the radial center has four coincident vertices).
 * 
 * \param cell Pyramid_Mesh (global) cell number
 *
 * \return vertices - the coordinates of the cell's eight vertices
 */
Pyramid_Mesh::vf_double Pyramid_Mesh::get_vertices(int cell) const
{
    Require (cell>0 && cell <= num_cells());
    Require (coord>get_dim() == 3);

    const int num_verts_face = 4;
    const int num_verts_cell = 8;
    const int loz_face = 5;
    const in hiz_face = 5;

    // get the vertices for the low z face of the cell
    vf_double cell_vertices = get_vertices(cell, loz_face);

    // get the vertices for the high z face of the cell
    vf_double hiz_face_vertices = get_vertices(cell,hiz_face);

    Require (cell_vertices.size() == coord->get_dim());
    Require (cell_vertices[0].size() == num_verts_face);
    Require (cell_vertices[1].size() == num_verts_face);
    Require (cell_vertices[2].size() == num_verts_face);

    Require (hiz_face_vertices.size() == coord->get_dim());
    Require (hiz_face_vertices[0].size() == num_verts_face);
    Require (hiz_face_vertices[1].size() == num_verts_face);
    Require (hiz_face_vertices[2].size() == num_verts_face);

    // add the hiz face vertices to the loz face vertices
    for (int v=0; v<coord->get_dim();v++)
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
 * \brief Calculate the vertices for a give face of a Pyramid_Mesh cell. 
 * 
 * During normal use of the Pyramid_Mesh, the explicit cell vertices are not
 * required. However, they are needed for graphics dumps.  The Pyramid_Mesh
 * is always in 3-D XYZ geometry and each face always has four vertices ( the
 * cell at the radial center has four coincident vertices).
 * 
 * \param cell Pyramid_Mesh (global) cell_number
 * \param face face of cell
 *
 * \return vertices of the coordinates of the four vertices defining a face
 */
Pyramid_mesh::vf_double Pyramid_Mesh::get_vertices(int cell, int face) const
{
    Require (face> 0 && face <= 6);
    Require (cell>0 && cell<= layout.num_cells());
    Require (coord->get_dim() == 3);

    vf_double face_vertices(coord->get_dim());
    sf_double single_vert(coord->get_dim(),0.0);

    double lox = get_low_x(cell);
    double hix = get_high_x(cell);
    double small_y = lox*tan_beta;
    double large_y = hix*tan_beta;
    double small_z = lox*tan_beta;
    double large_z = hix*tan_beta;

    // low x face or high x face
    if (face == 1 || face ==2)
    {
	double y_variable;
	double z_variable;

	if (face == 1)
	{
	    single_vert[0]=lox;
	    y_variable=small_y;
	    z_variable=small_z;
	}
	else if(face ==2)
	{
	    single_vert[0] = hix;
	    y_variable = large_y;
	    z_variable= large_z;
	}
	
	single_vert[1] = -y_variable;
	single_vert[2] = -z_variable;
	for(int i = 0; i <coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[1] = y_variable;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[2]= z_variable;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);


	single_vert[1]=-y_variable;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);
    }

    // low y face or high y face
    else if (face == 3 || face == 4)
    {
	double y_plusminus;
	if(face == 3)
	    y_plusminus = -1.0;
	else if(face == 4)
	    y_plusminus = 1.0;

	single_vert[0] =lox;
	single_vert[1]=y_plusminus*small_y;
	single_vert[2]=-small_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] =hix;
	single_vert[1]=y_plusminus*large_y;
	single_vert[2]=-large_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] =hix;
	single_vert[1]=y_plusminus*large_y;
	single_vert[2]=+large_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);	

	single_vert[0] =lox;
	single_vert[1]=y_plusminus*small_y;
	single_vert[2]=+small_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);	
    }

    // low z face or high z face
    
    else if (face == 5 || face == 6)
    {
	double z_plusminus;
	if(face == 5)
	    z_plusminus = -1.0;
	else if(face == 6)
	    z_plusminus = 1.0;

	single_vert[0] =lox;
	single_vert[1]=-small_y;
	single_vert[2]=z_plusminus*small_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] =hix;
	single_vert[1]=-large_y;
	single_vert[2]=z_plusminus*large_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);

	single_vert[0] =hix;
	single_vert[1]= large_y;
	single_vert[2]=z_plusminus*large_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);	

	single_vert[0] =lox;
	single_vert[1]=small_y;
	single_vert[2]=plus_minus*small_z;
	for(int i = 0; i<coord->get_dim(); i++)
	    face_vertices[i].push_back(single_vert[i]);	
    }

    // check that there are 3 dimensions a 4 vertices
    Ensure (face_vertices.size() == coord->get_dim());
    Ensure (face_vertices[0].size == 4);
    Ensure (face_vertices[1].size == 4);
    Ensure (face_vertices[2].size == 4);
 
    // return the four vertices for the face
    return face_vertices;
}
//---------------------------------------------------------------------------//
// Interface for graphics dumps
//---------------------------------------------------------------------------//
/*!
 * \brief Return turn cell type for each cell in the Pyramid_Mesh
 */
Pyramid_Mesh::sf_int Pyramid_Mesh::get_cell_types() const
{
    std::vector<int> cell_type(layout.num_cells());

    // all cells in a Pyramid_Mesh are general, 8-node hexahedrons
    std::fill(cell_type.begin(),cell_type.end(),
	      rtt_viz::eight_node_hexahedron);

    return cell_type;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the coordinates of all nodes in the mesh.
 * 
 * For each cell in the Pyramid_Mesh, the coordinates of all eight nodes are
 * returned. Thus, all interior nodes are replicated twice.  This approach,
 * although memory-inefficient, is easier to code especially considering that
 * the Pyramid_Mesh does not already have the raw vertex data.
 *
 */
Pyramid_Mesh::vf_double Pyramid_Mesh::get_point_coord() const
{
    using std::vector;

    // number of vertices is always 8; always 3D
    const in num_verts_cell = 8;
    int vert_index;
    Check (coord->get_dim() == 3);

    // weakly check the validity of num_cells()
    Check(num_cells() >0 );

    //initialize the return vector
    vector<vector<double> > return_coord(num_cells()*num_verts_cell);

    // for each cell, get_vertices and reverse the columns and rows
    for (int cell=1, cell<=num_cells(); cell++)
    {
	// get the cell's verices[dim][8]
	vector<vector<double>> cell_verts = Pyramid_Mesh::get_vertices(cell);

	// check validity of cell vertices vector
	Check (cell_verts.size() == coord->get_dim());

	// loop over all 8 nodes for this cell
	for (int node = 0; node <num_verts_cell; node++)
	{
	    //calculate running vertex index
	    vertex_index = (cell-1)*num_verts_cell+node;

	    //resize each vertice's return_coord to num dimensions
	    return_coord[vert_index].resize(coord->get_dim());

	    // re-assign point coordinates to return vector
	    for (int dim=0; dim<coord->get_dim(); dim++)
	    {
		Check (cell_verts[dim].size() == num_verts_cell);
		return_coord[vert_index][dim]=cell_verts[dim][node];
	    }
	}
    }

    // return the coordinate for all nodes of all cells
    return return_coord;
}





} // end namespace rtt_mc

#endif

//---------------------------------------------------------------------------//
//                 end of Pyramid_Mesh.cc
//---------------------------------------------------------------------------//
