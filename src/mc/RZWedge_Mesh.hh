//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/RZWedge_Mesh.hh
 * \author Todd J. Urbatsch
 * \date   Wed Apr  5 16:01:53 2000
 * \brief  RZWedge_Mesh header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_RZWedge_Mesh_hh__
#define __mc_RZWedge_Mesh_hh__

#include "Coord_sys.hh"
#include "AMR_Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "Math.hh"
#include "Constants.hh"
#include "Sampler.hh"
#include <vector>
#include <iostream>
#include <string>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class RZWedge_Mesh
 *
 * \brief An XYZ Wedge mesh constructed from an RZ mesh.

 * For Monte Carlo applications, we model an RZ mesh with a
 * three-dimensional, Cartesian wedge mesh where all cells are hexes except
 * for on the r=0 axis, where the inner face degenerates to a line.  The
 * x-direction is aligned with the r-direction.  The construction of the mesh
 * requires an RZ mesh and an unfolding angle, theta.  The radial dimensions
 * are adjusted for the particular theta in order to conserve volume.  The
 * mesh is one cell thick in the y-direction (or "azimuthal" direction) with
 * reflecting boundary conditions at the high and low y-sides (at +theta/2
 * and -theta/2).

 * The RZWedge_Mesh is fully AMR compatible and it uses the AMR_Layout class
 * internally to support this functionality.  The AMR_Layout must have the
 * following coarse face ordering: (lo x, hi x, lo y, hi y, lo z, hi z).
 * Across a coarse face the cells are ordered low to high.  Thus, across the
 * hi z coarse face index 1 must point to the low x cell and index 2 must
 * point to the high x cell.
 
 */
// revision history:
// -----------------
// 13 Apr 2000 : committed on Thursday, April 13, 2000.
//  5 Jul 2000 : added dependence on Coord_sys (needed in transport...)
// 17-JUL-2000 : added get_cell_pair function for graphics dumping.
// 24-JUL-2000 : replaced Layout with AMR_Layout
// 11-AUG-2000 : fixed error in sample_pos() (tilt), the z position was
//               incorrectly being weighted by the z dimension; also added
//               in_cell() test function to determine if a set of coordinates 
//               is in a cell
// 
//===========================================================================//

class RZWedge_Mesh 
{
  public:
    // typedefs used throughout RZWedge_Mesh class
    typedef rtt_dsxx::SP<RZWedge_Mesh>        SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>           SP_Coord;
    typedef rtt_rng::Sprng                    rng_Sprng;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::string                       std_string;

    // forward declarations of cell-centered fields
    template<class T> class CCSF;
    template<class T> class CCVF;   

    // Handy typedefs to CC fields (not formally needed in KCC3.3+).
    typedef CCSF<double>     CCSF_double;
    typedef CCSF<int>        CCSF_int;
    typedef CCVF<double>     CCVF_double;
    typedef CCVF<int>        CCVF_int;
    typedef CCSF<std_string> CCSF_string;

  private:
    // base class reference to a derived coordinate system class
    // (the RZWedge_Mesh is always three-dimensional, Cartesian)
    SP_Coord coord;

    // layout (cell connectivity) of mesh
    AMR_Layout layout;

    // vector<vector> of x- and z-extents of each cell
    vf_double cell_xz_extents;

    // Indicator whether this is a submesh 
    bool submesh;

    // wedge angle data; precalculated 
    double theta_radians;
    double theta_degrees;
    double tan_half_theta;
    double sin_half_theta;
    double cos_half_theta;

    // >>> private implementations

    // function to calculate frequently used wedge data
    void calc_wedge_angle_data(const double);
   
    // Private copy and assignment operators (can't copy or assign a mesh).
    RZWedge_Mesh(const RZWedge_Mesh &);
    RZWedge_Mesh& operator=(const RZWedge_Mesh &);

  public:
    // Constructor.
    RZWedge_Mesh(SP_Coord, AMR_Layout &, vf_double &, double, bool = false); 

    // >>> Member functions used by RZWedge_Mesh-dependent classes

    //! Return the number of cells.
    int num_cells() const { return layout.num_cells(); }

    // Return the x- and z-dimension cell extents.
    double get_low_x(int cell)  const { return cell_xz_extents[cell-1][0]; }
    double get_high_x(int cell) const { return cell_xz_extents[cell-1][1]; }
    double get_low_z(int cell)  const { return cell_xz_extents[cell-1][2]; }
    double get_high_z(int cell) const { return cell_xz_extents[cell-1][3]; }

    // Get the midpoint of a cell for a given dimension.
    inline double get_x_midpoint(int cell) const;
    inline double get_y_midpoint(int cell) const;
    inline double get_z_midpoint(int cell) const;

    // get the dimension of a cell for a given coordinate 
    // (y-coord => dim at midpoint)  
    inline double dim(const int coordinate, const int cell) const;

    // Determine if a position is in a cell.
    bool in_cell(int, const sf_double &) const;

    // Diagnostic functions.
    void print(std::ostream &) const;
    void print(std::ostream &, int) const;

    // References to embedded objects.
    const AMR_Layout& get_Layout()  const { return layout; }
    const Coord_sys&  get_Coord()   const { return *coord; }
    SP_Coord          get_SPCoord() const { return coord; }

    // Services required for graphics dumps.
    sf_int get_cell_types() const;
    vf_double get_point_coord() const;
    vf_int get_cell_pair() const;

    // Required services for transport and source.
    int next_cell(int cell, int face) const { return layout(cell, face, 1); }
    int get_cell(const sf_double &) const;
    double get_db(const sf_double &, const sf_double &, int, int &) const;
    inline sf_double get_normal(int, int) const;
    inline sf_double get_normal_in(int, int) const;
    inline double volume(int) const;
    inline double face_area(int, int) const;
    vf_double get_vertices(int) const;
    vf_double get_vertices(int, int) const;
    inline sf_double sample_pos(int, rng_Sprng &) const;
    inline sf_double sample_pos(int, rng_Sprng &, sf_double, double) const;
    inline sf_double sample_pos_on_face(int, int, rng_Sprng &) const;
    int get_bndface(std_string, int) const;
    sf_int get_surcells(std_string) const;
    bool check_defined_surcells(const std_string, const sf_int &) const;
    inline sf_int get_neighbors(int) const;
    bool full_Mesh() const { return !submesh; }

    // Overloaded operators.
    bool operator==(const RZWedge_Mesh &) const;
    bool operator!=(const RZWedge_Mesh &rhs) const { return !(*this == rhs); }

};

//---------------------------------------------------------------------------//
// RZWEDGE_MESH INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream &output, const RZWedge_Mesh &object);

//---------------------------------------------------------------------------//
// RZWEDGE_MESH SERVICES FOR RZWEDGE_MESH-DEPENDENT CLASSES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the midpoint of a cell for the x-dimension.
 *
 * \param cell RZWedge_Mesh cell.
 *
 * \return midpoint of x-dimension
 */
double RZWedge_Mesh::get_x_midpoint(int cell) const
{
    Check ( (cell > 0) && (cell <= num_cells()) );

    return 0.5 * (get_low_x(cell) + get_high_x(cell));
}
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the midpoint of a cell for the y-dimension.
 *
 * \param cell RZWedge_Mesh cell.
 *
 * \return midpoint of y-dimension
 */
double RZWedge_Mesh::get_y_midpoint(int cell) const
{ 
    Check ( (cell > 0) && (cell <= num_cells()) );

    return 0.0; 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the midpoint of a cell for the z-dimension.
 *
 * \param cell RZWedge_Mesh cell.
 *
 * \return midpoint of z-dimension
 */
double RZWedge_Mesh::get_z_midpoint(int cell) const
{
    Check ( (cell > 0) && (cell <= num_cells()) );

    return 0.5 * (get_low_z(cell) + get_high_z(cell));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the dimension of the cell in the requested coordinate.
 *
 * \param cell  RZWedge_Mesh cell.
 * \param coord coord
 *
 * \return dim the dimension (length) of the cell in the requested coordinate
 */
double RZWedge_Mesh::dim(const int coordinate, const int cell) const
{
    Check ( (cell > 0) && (cell <= num_cells()) );
    Check ( (coordinate > 0) && (coordinate <= 3) );

    // return value
    double dimension = 0.0;

    // x-coordinate
    if (coordinate == 1)
	dimension = get_high_x(cell) - get_low_x(cell);

    // y-coordinate
    else if (coordinate == 2)
	dimension = 2.0 * get_x_midpoint(cell) * tan_half_theta;

    // z-coordinate
    else if (coordinate == 3)
	dimension = get_high_z(cell) - get_low_z(cell); 
    
    else 
	Insist (0, "Requested coordinate in rzwedge_mesh's dim not valid!"); 
    
    return dimension;
}

//---------------------------------------------------------------------------//
// RZWEDGE_MESH GENERALIZED MT SERVICES REQUIRED BY IMC
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the outward normal for the particular face of a cell.
 *
 * \param cell cell -- not used -- each RZWedge_Mesh cell has same normals.
 * \param face face of cell for which to return outward normal.
 *
 * \return outward normal
 */
RZWedge_Mesh::sf_double RZWedge_Mesh::get_normal(int cell, int face) const
{
    Check (coord->get_dim() == 3);
    Check ((face >= 1) && (face <= 6));

    sf_double normal(coord->get_dim(), 0.0);

    // low x face
    if (face == 1)
	normal[0] = -1.0;

    // high x face
    else if (face == 2)
	normal[0] = 1.0;

    // low y face
    else if (face == 3)
    {
	normal[0] = -sin_half_theta;
	normal[1] = -cos_half_theta;
    }

    // high y face
    else if (face == 4)
    {
	normal[0] = -sin_half_theta;
	normal[1] =  cos_half_theta;
    }

    // low z face
    else if (face == 5)
	normal[2] = -1.0;

    // high z face
    else if (face == 6)
	normal[2] = 1.0;

    // return outward normal
    return normal;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the inward normal for the particular face of a cell.
 *
 * \param cell cell -- not used -- each RZWedge_Mesh cell has same normals.
 * \param face face of cell for which to return inward normal.
 *
 * \return inward normal
 */
RZWedge_Mesh::sf_double RZWedge_Mesh::get_normal_in(int cell, int face) const
{
    Check (coord->get_dim() == 3);
    Check ((face >= 1) && (face <= 6));

    // initialize inward normal
    sf_double normal_in(coord->get_dim(), 0.0);

    // get outward normal first
    sf_double normal = get_normal(cell, face);
    Check (normal.size() == coord->get_dim());

    // reverse direction
    for (int dir = 0; dir < coord->get_dim(); dir++)
	normal_in[dir] = -normal[dir];

    // return inward normal
    return normal_in;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the volume of an RZWedge_Mesh cell.
 *
 * \param cell RZWedge_Mesh cell.
 *
 * RZWedge cell construction: align R with X; unfold RZ mesh about X;
 * correct X values to conserve radial/XY area.
 *
 * \return volume of wedge cell
 */
double RZWedge_Mesh::volume(int cell) const
{
    Require (cell > 0 && cell <= num_cells());

    double lox = get_low_x(cell);
    double hix = get_high_x(cell);
    
    double vol = tan_half_theta * (hix - lox) * (hix + lox) *
	(get_high_z(cell) - get_low_z(cell));

    Ensure (vol > 0.0);

    return vol;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the area of a face in an RZWedge_Mesh cell.
 *
 * \param cell RZWedge_Mesh cell.
 * \param face face of cell.
 *
 * RZWedge cell construction: align R with X; unfold RZ mesh about X;
 * correct X values to conserve radial/XY area.
 *
 * \return volume of wedge cell
 */
double RZWedge_Mesh::face_area(int cell, int face) const
{
    Require (face > 0 && face <= 6);
    Require (cell > 0 && cell <= num_cells());

    // initialized face area
    double face_area = 0.0;

    // low x face
    if (face == 1)
	face_area = 2.0 * get_low_x(cell) * tan_half_theta * 
	    (get_high_z(cell) - get_low_z(cell));

    // high x face
    else if (face == 2)
	face_area = 2.0 * get_high_x(cell) * tan_half_theta * 
	    (get_high_z(cell) - get_low_z(cell));

    // low y face or high y face
    else if (face == 3 || face == 4)
	face_area = (get_high_x(cell) - get_low_x(cell)) /
	    cos_half_theta * (get_high_z(cell) - get_low_z(cell));

    // low z face or high z face 
    else if (face == 5 || face ==6)
	face_area = (get_high_x(cell) - get_low_x(cell)) *
	    (get_high_x(cell) + get_low_x(cell)) * tan_half_theta;

    // return face area 
    return face_area;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a position uniformly in an RZWedge_Mesh cell.
 *
 * \param cell RZWedge_Mesh cell.
 * \param rng_Sprng random number object
 *
 * \return location in cell.
 */
RZWedge_Mesh::sf_double RZWedge_Mesh::sample_pos(int cell, 
						 rng_Sprng &random) const
{
    Require (cell > 0 && cell <= num_cells());
    Check (coord->get_dim() == 3);

    // initialize location vector
    sf_double position(coord->get_dim(), 0.0);

    // get x-dimension cell extents
    double lox =  get_low_x(cell);
    double hix = get_high_x(cell);

    // calculate y-values corresponding to x-dimension cell extents
    double loy = lox * tan_half_theta;
    double hiy = hix * tan_half_theta;

    // sample x uniformly in wedge area of x-y face
    position[0] = rtt_mc::sampler::sample_general_linear(random, lox, hix,
							 loy, hiy);

    // calculate y-values corresponding to sampled x-value
    double pos_y = position[0] * tan_half_theta;
    double neg_y = -pos_y;
    Check (pos_y >= loy && pos_y <= hiy);

    // sample y uniformly
    position[1] = neg_y + random.ran() * (pos_y - neg_y);

    // get z-dimension cell extents
    double loz = get_low_z(cell);
    double hiz = get_high_z(cell);

    // sample z uniformly
    position[2] = loz + random.ran() * (hiz - loz);

    // Check that sampled position is within cell
    Check (position[0] >= lox   && position[0] <= hix);
    Check (position[1] >= neg_y && position[1] <= pos_y);
    Check (position[2] >= loz   && position[2] <= hiz);

    // return location vector
    return position;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a position in an RZWedge_Mesh cell according to a slope.
 *
 * pdf is linear and assumed separable in each dimension: 
 * f(x,y,x) = f(x)f(y)f(z).
 *
 * \param cell RZWedge_Mesh cell.
 * \param rng_Sprng random number object
 * \param slope slope of sampling weight function
 * \param center_pt cell-centered value of sampling weight function
 *
 * \return location in cell.
 */
RZWedge_Mesh::sf_double RZWedge_Mesh::sample_pos(int cell, 
						 rng_Sprng &random,
						 sf_double slope, 
						 double center_pt ) const
{
    using global::max;

    Require (cell > 0 && cell <= num_cells());
    Check (coord->get_dim() == 3);

    // initialize location vector
    sf_double position(coord->get_dim(), 0.0);

    // get x-dimension cell extents
    double lox       = get_low_x(cell);
    double hix       = get_high_x(cell);
    double half_delx = 0.5 * (hix - lox);

    // calculate weighting function values at x-dimension cell extents
    double half_delw = slope[0]  * half_delx;
    double low_w     = center_pt - half_delw;
    double high_w    = center_pt + half_delw;

    // set low_w and high_w to zero if it is less than zero
    if (low_w < 0)
    {
	Check (low_w > -1.e-12);
	low_w = 0.0;
    }
    if (high_w < 0)
    {
	Check (high_w > -1.e-12);
	high_w = 0.0;
    }

    // calculate posititve y-values corresponding to x-dimension cell extents
    double loy = lox * tan_half_theta;
    double hiy = hix * tan_half_theta;

    // calc extents of function (wt-fn * half_y_value) from which to sample
    double lof = loy * low_w;
    double hif = hiy * high_w;

    // sample x 
    position[0] = rtt_mc::sampler::sample_general_linear(random, lox, hix,
							 lof, hif);

    // calculate y-values corresponding to sampled x-value
    double pos_y = position[0] * tan_half_theta;
    double neg_y = -pos_y;
    Check (pos_y >= loy && pos_y <= hiy);

    // sample y uniformly
    Check (rtt_mc::global::soft_equiv(slope[1],0.0)); 
    position[1] = neg_y + random.ran() * (pos_y - neg_y);

    // get z-dimension cell extents
    double loz       = get_low_z(cell);
    double hiz       = get_high_z(cell);
    double half_delz = 0.5 * (hiz - loz);

    // calculate weighting function at z-dimension cell extents
    half_delw = slope[2] * half_delz;
    low_w     = center_pt - half_delw;
    high_w    = center_pt + half_delw;

    // set low_w and high_w to zero if it is less than zero
    if (low_w < 0)
    {
	Check (low_w > -1.e-12);
	low_w = 0.0;
    }
    if (high_w < 0)
    {
	Check (high_w > -1.e-12);
	high_w = 0.0;
    }

    // calculate extents of function from which to sample
    lof = low_w;
    hif = high_w;

    // sample z 
    position[2] = rtt_mc::sampler::sample_general_linear(random, loz, hiz,
							 lof, hif);

    // Check that sampled position is within cell
    Check (position[0] >= lox   && position[0] <= hix);
    Check (position[1] >= neg_y && position[1] <= pos_y);
    Check (position[2] >= loz   && position[2] <= hiz);

    // return location vector
    return position;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a position on an RZWedge_Mesh cell face.
 *
 * \param cell RZWedge_Mesh cell.
 * \param face RZWedge_Mesh face.
 * \param rng_Sprng random number object
 *
 * \return location in cell.
 */
RZWedge_Mesh::sf_double RZWedge_Mesh::sample_pos_on_face(int cell, int face,
							 rng_Sprng &random) const
{
    Require (cell > 0 && cell <= num_cells());
    Require (face > 0 && face <= 6);
    Check (coord->get_dim() == 3);

    // initialize position vector
    sf_double position(coord->get_dim(), 0.0);

    // high x face
    if (face == 2)
    {
	double hix = get_high_x(cell);
	double hiy = tan_half_theta * hix;
	double loz = get_low_z(cell);
	double hiz = get_high_z(cell);

	position[0] = hix;
	position[1] = -hiy + (random.ran() * 2.0 * hiy);
	position[2] = loz + random.ran() * (hiz - loz);
    }

    // low z face  or  high z face
    else if (face == 5 || face == 6)
    {
	double lox = get_low_x(cell);
	double hix = get_high_x(cell);
	double loy = tan_half_theta * lox;
	double hiy = tan_half_theta * hix;

	position[0] = rtt_mc::sampler::sample_general_linear(random, lox,
							     hix, loy, hiy);
	double y = tan_half_theta * position[0];
	position[1] = -y + (random.ran() * 2.0 * y);

	if (face == 5)
	    position[2] = get_low_z(cell);
	else
	    position[2] = get_high_z(cell);
    }

    else
    {
	position[0] = 0.0;
	Insist (0, "Should not be sampling lox/loy/hiy RZWedge_Mesh face!"); 
    }

    // return the sampled position
    return position;
}
//---------------------------------------------------------------------------//
// return a sf_int list of a cells neighbors

RZWedge_Mesh::sf_int RZWedge_Mesh::get_neighbors(int cell) const
{
    sf_int neighbors(layout.num_faces(cell));
    
    for (int face = 1; face <= neighbors.size(); face++)
	neighbors[face-1] = layout(cell, face, 1);

    return neighbors;
}

//===========================================================================//
// class RZWedge_Mesh::CCSF (copied from class OS_Mesh::CCSF)
//
// cell-centered scalar fields
// Note: we can't build empty fields (ie. the mesh has to have at least 1
// cell 
//===========================================================================//

template<class T>
class RZWedge_Mesh::CCSF
{
  public:
    // STL style typedefs.
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef typename std::vector<T>::pointer pointer;
    typedef typename std::vector<T>::const_pointer const_pointer;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::size_type size_type;

  private:
    // SP back to RZWedge_Mesh.
    SP_Mesh mesh;

    // Data in field, (num_cells).
    std::vector<T> data;

  public:
    // Explicit constructor.
    inline explicit CCSF(SP_Mesh);

    // Additional constructors.
    inline CCSF(SP_Mesh, const std::vector<T> &);

    // Return reference to mesh.
    const RZWedge_Mesh& get_Mesh() const { return *mesh; }

    // Subscripting.

    // We use () to indicate absolute cell number and [] to indicate
    // vector-style (0:n-1) indexing.
    const_reference operator()(int cell) const { return data[cell-1]; }
    reference operator()(int cell) { return data[cell-1]; }
 
    const_reference operator[](int index) const { return data[index]; }
    reference operator[](int index) { return data[index]; }

    // STL style functions.
    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }
    
    iterator end() { return data.end(); }
    const_iterator end() const {return data.end();}
    
    size_type size() const { return data.size(); }
    bool empty() const { return data.empty(); }
}; 

//---------------------------------------------------------------------------//
// RZWedge_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//
// CCSF explicit constructor

template<class T>
RZWedge_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_) 
    : mesh(mesh_), data(mesh->num_cells()) 
{
    Require (mesh);
    Ensure  (!empty());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
RZWedge_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_, const std::vector<T> &array)
    : mesh(mesh_), data(array)
{
    Require (mesh)
    Ensure  (data.size() == mesh->num_cells());
    Ensure  (!empty());
}

//===========================================================================//
// class RZWedge_Mesh::CCVF (copied from class OS_Mesh::CCVF)
//
// cell-centered vector fields
//===========================================================================//

template<class T>
class RZWedge_Mesh::CCVF
{
  public:
    // STL style typedefs.
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef typename std::vector<T>::pointer pointer;
    typedef typename std::vector<T>::const_pointer const_pointer;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::size_type size_type;
    
  private:
    // SP back to RZWedge_Mesh
    SP_Mesh mesh;
    
    // 2-D field vector, (dimension, num_cells).
    std::vector<std::vector<T> > data;

  public:
    // Inline explicit constructor.
    inline explicit CCVF(SP_Mesh);

    // Additional constructors.
    inline CCVF(SP_Mesh, const std::vector<std::vector<T> > &);

    // Return reference to mesh.
    const RZWedge_Mesh& get_Mesh() const { return *mesh; }

    // Subscripting.
    inline const T& operator()(int, int) const;
    inline T& operator()(int, int);

    // Getting a CC vector in a cell.
    inline std::vector<T> operator()(int) const;

    // STL style functions.
    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }
    inline iterator begin(int i);
    inline const_iterator begin(int i) const;
    
    iterator end() { return data.end(); }
    const_iterator end() const { return data.end();}
    inline iterator end(int i);
    inline const_iterator end(int i) const;
    
    size_type size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    inline size_type size(int i) const;
    inline bool empty(int i) const;
}; 

//---------------------------------------------------------------------------//
// RZWedge_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//
// CCVF explicit constructor

template<class T>
RZWedge_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_)
    : mesh(mesh_), data(mesh->get_Coord().get_dim())
{
    Require (mesh);
    Check   (mesh->get_Coord().get_dim() == 3);

    // initialize data array
    for (int i = 0; i < mesh->get_Coord().get_dim(); i++)
	data[i].resize(mesh->num_cells());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
RZWedge_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_, 
			    const std::vector<std::vector<T> > &array)
    : mesh(mesh_), data(array)
{
    // check things out
    Ensure (data.size() == mesh->get_Coord().get_dim());
    Ensure (mesh->get_Coord().get_dim() == 3);

    for (int dim = 0; dim < mesh->get_Coord().get_dim(); dim++)
	Ensure (data[dim].size() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
const T& RZWedge_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim-1][cell-1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
T& RZWedge_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim-1][cell-1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
std::vector<T> RZWedge_Mesh::CCVF<T>::operator()(int cell) const
{
    // declare return vector
    std::vector<T> x;
    
    // loop through dimensions and make return vector for this cell
    for (int i = 0; i < data.size(); i++)
	x.push_back(data[i][cell-1]);

    // return
    Ensure (x.size() == data.size());
    return x;
} 

//---------------------------------------------------------------------------//
// STL style functionality for CCVF fields

template<class T>
RZWedge_Mesh::CCVF<T>::iterator RZWedge_Mesh::CCVF<T>::begin(int i)
{
    Require(i > 0 && i <= data.size());
    return data[i-1].begin();
}

template<class T>
RZWedge_Mesh::CCVF<T>::const_iterator RZWedge_Mesh::CCVF<T>::begin(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].begin();
}

template<class T>
RZWedge_Mesh::CCVF<T>::iterator RZWedge_Mesh::CCVF<T>::end(int i)
{
    Require(i > 0 && i <= data.size());
    return data[i-1].end();
}

template<class T>
RZWedge_Mesh::CCVF<T>::const_iterator RZWedge_Mesh::CCVF<T>::end(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].end();
}

template<class T>
RZWedge_Mesh::CCVF<T>::size_type RZWedge_Mesh::CCVF<T>::size(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].size();
}

template<class T>
bool RZWedge_Mesh::CCVF<T>::empty(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].empty();
}

} // end namespace rtt_mc

#endif                          // __mc_RZWedge_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of mc/RZWedge_Mesh.hh
//---------------------------------------------------------------------------//
