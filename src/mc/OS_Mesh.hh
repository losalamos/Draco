//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/OS_Mesh.hh
 * \author Thomas M. Evans
 * \date   Tue Feb  3 16:50:13 1998
 * \brief  OS_Mesh class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_OS_Mesh_hh__
#define __mc_OS_Mesh_hh__

#include "Coord_sys.hh"
#include "Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class OS_Mesh
 *
 * \brief Orthogonal-Structured Mesh Class for Monte Carlo applications.
 *
 * The OS_Mesh class is a generalized orthogonal-structured mesh class for
 * Monte Carlo applications.
 */
/*!
 * \example mc/test/tstOSMesh.cc
 *
 * Example usage of the OS_Mesh class services for Monte Carlo transport.
 */
// revision history:
// -----------------
//  0) original
//  1)   2-6-98 : changed storage of Layout to a full object instead of a 
//                smart pointer
//  2)  2-20-98 : added assertion to constructor to make sure that the number
//                of cells defined in the mesh is consistent with the Layout
//  3)  3-16-98 : added Get_normal public member function, OS_Meshes know the
//                normals a priori so they can just be sent.  Non_OS meshes
//                may need a yet undefined (Calc_normal) function from
//                Coord_sys to calculate the normal based on the cell
//                vertices
//  4)  3-18-98 : added generalized Mesh constructor based on cell vertices,
//                use this in the Get_db and other member functions
//  5)  4-14-98 : added get_vertex and get_cellpair to compliment get_Coord 
//                and get_Layout; these functions are needed in
//                Parallel_Builder 
//  6)  5-14-98 : added overloaded operator== for Checks and design by
//                contract purposes
//  7)  6-10-98 : added constructors to CCSF and CCVF for taking a vector to
//                initialize the field
//  8)  6-11-98 : added get_normal_in public member function; slightly 
//                modified from get_normal to give the inward normal
//  9)  6-12-98 : added piecewise-linear sampling to sample_pos member
//                function
// 10)  6-15-98 : added vector return operator() to CCVF for getting a CC
//                vector from a cell
// 11)   7-6-98 : added submesh indicator for domain decomposition problems
// 12)   7-8-98 : added new get_vertices(int cell) member   
// 13)  4-13-99 : moved to mc package
// 14)   7-6-99 : moved CC field definitions to a location after the OS_Mesh
//                definition to improve readability; made the CCSF class STL
//                compliant (ie. it now has iterators that can be used to 
//                navigate through the class)-->we can now use the algorithms 
//                provided by the STL to do operations on CCSF fields, also
//                this move presupposes the use of Randy's matprops class
// 15) 18-NOV-99: added get_neighbors member function that returns a
//                vector<int> that contains the neighbors of a cell.
// 16) 29-NOV-99: added full_Mesh() function that returns the value of
//                submesh, used for determining whether a mesh is a master
//                (global) mesh or some sort of decomposed 
// 17) 20-DEC-99: added STL random access iterator functionality to CCVF
// 18) 27-JAN-00: added get_cell_types for graphics data dumping
// 19) 18-AUG-00: added new signature for next_cell function
// 20) 31-AUG-00: added get_spatial_dimension() function
// 
//===========================================================================//
    
class OS_Mesh
{
  public:
    // Useful typdefs to std:: namespace members.
    typedef rtt_dsxx::SP<OS_Mesh>             SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>           SP_Coord_sys;
    typedef rtt_rng::Sprng                    rng_Sprng;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::string                       std_string;

    // Forward Declarations of cell centered fields.
    template<class T> class CCSF;
    template<class T> class CCVF;
   
    // Handy typedefs to CC fields (not formally needed in KCC3.3+).
    typedef CCSF<double>     CCSF_double;
    typedef CCSF<int>        CCSF_int;
    typedef CCVF<double>     CCVF_double;
    typedef CCVF<int>        CCVF_int;
    typedef CCSF<std_string> CCSF_string;

  private:
    // Base class reference to a derived coord class.
    rtt_dsxx::SP<Coord_sys> coord;

    // Layout of mesh.
    Layout layout;

    // Vertices in mesh.
    vf_double vertex;

    // Cell-pairings of cell to its vertices.
    vf_int cell_pair;

    // Area of surfaces on each dimension.
    vf_double sur;

    // Indicator whether this is a submesh.
    bool submesh;

    // >>> Private implementations.

    // Calculate a surface array from the vertices of the mesh.
    void calc_surface();

    // Private copy and assignment operators; can't copy or assign a mesh.
    OS_Mesh(const OS_Mesh &);
    OS_Mesh& operator=(const OS_Mesh &);

  public:
    // Generalized constructor for all mesh types.
    OS_Mesh(SP_Coord_sys, Layout &, vf_double &, vf_int &, bool = false);

    // >>> Member functions used by the OS_Mesh-dependent classes.

    // Mesh dimensionality functions.

    // Give the dimension and begin and end return the beginning and ending
    // coordinate along that dimension.
    inline double begin(int) const;
    inline double end(int) const;

    // Return number of cells.
    int num_cells() const { return layout.num_cells(); }

    // Cell dimensionality functions.

    // Find minimum and maximum dimension of cell.
    inline double min(int, int) const;
    inline double max(int, int) const;

    // Find centerpoint of cell and width of cell.
    inline double pos(int, int) const;
    double dim(int d, int cell) const { return max(d, cell) - min(d, cell); } 

    // Diagnostic functions.
    void print(std::ostream &) const;
    void print(std::ostream &, int) const;

    // References to imbedded objects and data required for
    // Parallel_Building. 
    const Layout& get_Layout() const { return layout; }
    const Coord_sys& get_Coord() const { return *coord; }
    rtt_dsxx::SP<Coord_sys> get_SPCoord() const { return coord; }
    const vf_double& get_vertex() const { return vertex; }
    const vf_int& get_cell_pair() const { return cell_pair; }

    // >>> Services required by ALL mesh types used in IMC-codes.

    // Services required for Graphics dumps.
    sf_int get_cell_types() const;
    vf_double get_point_coord() const;

    // Required services for transport and source.
    int get_spatial_dimension() const { return coord->get_dim(); }
    inline int next_cell(int, int, const sf_double & = sf_double()) const;
    int get_cell(const sf_double &) const;
    double get_db(const sf_double &, const sf_double &, int, int &) const; 
    inline sf_double get_normal(int, int) const;
    inline sf_double get_normal_in(int, int) const;
    inline double volume(int) const;
    inline double face_area(int, int) const;
    sf_int get_surcells(std_string) const;
    bool check_defined_surcells(const std_string, const sf_int &) const;
    int get_bndface(std_string, int) const;
    inline vf_double get_vertices(int, int) const;
    inline vf_double get_vertices(int) const;
    inline sf_double sample_pos(int, rng_Sprng &) const;
    inline sf_double sample_pos(int, rng_Sprng &, sf_double, double) const; 
    inline sf_double sample_pos_on_face(int, int, rng_Sprng &)	const;
    inline sf_int get_neighbors(int) const;
    bool full_Mesh() const { return !submesh; }

    // Overloaded operators.
    bool operator==(const OS_Mesh &) const;
    bool operator!=(const OS_Mesh &rhs) const { return !(*this == rhs); }
};

//===========================================================================//
// OS_MESH INLINE FUNCTIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream &output, const OS_Mesh &object);

//---------------------------------------------------------------------------//
// OS_MESH SERVICES FOR OS_MESH DEPENDENT CLASSES
//---------------------------------------------------------------------------//

double OS_Mesh::begin(int d) const 
{
    // find the minimum surface for d over the whole mesh
    return *std::min_element(vertex[d-1].begin(), vertex[d-1].end()); 
}

//---------------------------------------------------------------------------//

double OS_Mesh::end(int d) const 
{
    // find the maximum surface for d over the whole mesh
    return *std::max_element(vertex[d-1].begin(), vertex[d-1].end()); 
}

//---------------------------------------------------------------------------//
// find center position of cell

double OS_Mesh::pos(int d, int cell) const
{
    // set return value
    double return_pos = 0.0;
	
    // loop over all vertices and take average value to get the center
    // point 
    for (int i = 0; i < cell_pair[cell-1].size(); i++)
	return_pos += vertex[d-1][cell_pair[cell-1][i]-1];

    // return value
    return return_pos / static_cast<double>(cell_pair[cell-1].size());     
}

//---------------------------------------------------------------------------//
// find minimum dimension along d of cell

double OS_Mesh::min(int d, int cell) const 
{	
    // loop over all vertices and find the minimum
    double minimum = vertex[d-1][cell_pair[cell-1][0]-1];
    for (int i = 1; i < cell_pair[cell-1].size(); i++)
    {
	double point = vertex[d-1][cell_pair[cell-1][i]-1];

	// update the minimum value point
	if (point < minimum)
	    minimum = point;
    }
	
    // return minimum dimension
    return minimum;
}

//---------------------------------------------------------------------------//
// find maximum dimension of cell

double OS_Mesh::max(int d, int cell) const
{
    // loop over all vertices and find the maximum
    double maximum = vertex[d-1][cell_pair[cell-1][0]-1];
    for (int i = 1; i < cell_pair[cell-1].size(); i++)
    {
	double point = vertex[d-1][cell_pair[cell-1][i]-1];

	// update the maximum value point
	if (point > maximum)
	    maximum = point;
    }

    // return maximum dimension
    return maximum;
}

//---------------------------------------------------------------------------//
// OS_MESH GENERALIZED MT SERVICES REQUIRED BY IMC
//---------------------------------------------------------------------------//
// Determine the next cell, since OS_Mesh is a zero level mesh, the next cell 
// function will always neglect the r argument

int OS_Mesh::next_cell(int cell, int face, const sf_double &r) const
{
    Require (cell > 0 && cell <= layout.num_cells());
    return layout(cell, face);
}

//---------------------------------------------------------------------------//
// calculate volume of cell

double OS_Mesh::volume(int cell) const 
{
    // loop through dimensions and get volume
    double volume = 1.0;
    for (int d = 1; d <= coord->get_dim(); d++)
	volume *= dim(d, cell);

    // return volume
    return volume;
}

//---------------------------------------------------------------------------//
// calculate area of face on cell

double OS_Mesh::face_area(int cell, int face) const 
{
    // loop through dimensions and multiply off-dimension widths
    double face_area = 1.0;
    int dim_face_on = (face + 1) / 2;
    for (int d = 1; d <= coord->get_dim(); d++)
    {
	if (d != dim_face_on)
	    face_area *= dim(d, cell);
    }

    // return face_area
    return face_area;
}

//---------------------------------------------------------------------------//

OS_Mesh::sf_double OS_Mesh::get_normal(int cell, int face) const
{
    // OS_Meshes do not require any functionality from Coord_sys to 
    // calculate the outward normal, do simple return

    // normal always has 3 components, use Get_sdim()
    sf_double normal(coord->get_sdim(), 0.0);
	
    // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
    // one coordinate is non-zero    
    normal[(face-1)/2] = pow(-1.0, face);

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//

OS_Mesh::sf_double OS_Mesh::get_normal_in(int cell, int face) const
{
    // OS_Meshes do not require any functionality from Coord_sys to 
    // calculate the inward normal, do simple return

    // normal always has 3 components, use Get_sdim()
    sf_double normal(coord->get_sdim(), 0.0);
	
    // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
    // one coordinate is non-zero    
    normal[(face-1)/2] = pow(-1.0, face-1);

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell face

OS_Mesh::vf_double OS_Mesh::get_vertices(int cell, int face) const
{
    // determine the vertices along a cell-face

    // return vertices
    vf_double ret_vert(coord->get_dim());

    // determine axis dimension of surface (x=1, y=2, z=3)
    int axis = (face + 1)/2;
    double plane;
    if (2*axis - 1 == face)
	plane = min(axis, cell);
    else
	plane = max(axis, cell);

    // loop over vertices in cell and get the vertices that are in the plane
    for (int i = 0; i < cell_pair[cell-1].size(); i++)
	if (plane == vertex[axis-1][cell_pair[cell-1][i]-1])
	    for (int d = 0; d < coord->get_dim(); d++)
		ret_vert[d].push_back(vertex[d][cell_pair[cell-1][i]-1]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell

OS_Mesh::vf_double OS_Mesh::get_vertices(int cell) const
{
    // determine the vertices bounding a cell
    
    // return vertices
    vf_double ret_vert(coord->get_dim());

    // loop over cell vertices and build the cell vertices
    for (int i = 0; i < cell_pair[cell-1].size(); i++)
	for (int d = 0; d < coord->get_dim(); d++)
	    ret_vert[d].push_back(vertex[d][cell_pair[cell-1][i]-1]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// sample the position uniformly in a cell

OS_Mesh::sf_double OS_Mesh::sample_pos(int cell, rng_Sprng &random) const
{
    // assign minimums and maximums for cell dimensions
    sf_double vmin(coord->get_dim());
    sf_double vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    sf_double r = coord->sample_pos(vmin, vmax, random);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample the position in a cell given a tilt (or other slope-like function)

OS_Mesh::sf_double OS_Mesh::sample_pos(int cell, rng_Sprng &random,
				       sf_double slope, 
				       double center_pt) const
{
    // assign minimums and maximums for cells dimensions
    sf_double vmin(coord->get_dim());
    sf_double vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    sf_double r = coord->sample_pos(vmin, vmax, random, slope, center_pt);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample a position on a face

OS_Mesh::sf_double OS_Mesh::sample_pos_on_face(int cell, int face, 
					       rng_Sprng &random) const
{
    // assign minimums and maximums for cell dimensions
    sf_double vmin(coord->get_dim());
    sf_double vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    sf_double r = coord->sample_pos_on_face(vmin, vmax, face, random);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// return a sf_int list of a cells neighbors

OS_Mesh::sf_int OS_Mesh::get_neighbors(int cell) const
{
    sf_int neighbors(layout.num_faces(cell));
    
    for (int face = 1; face <= neighbors.size(); face++)
	neighbors[face-1] = layout(cell, face);

    return neighbors;
}

//===========================================================================//
// class OS_Mesh::CCSF
//
// cell-centered scalar fields
// Note: we can't build empty fields (ie. the mesh has to have at least 1
// cell 
//===========================================================================//

template<class T>
class OS_Mesh::CCSF
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
    // SP back to OS_Mesh.
    SP_Mesh mesh;

    // Data in field, (num_cells).
    std::vector<T> data;

  public:
    // Explicit constructor.
    inline explicit CCSF(SP_Mesh);

    // Additional constructors.
    inline CCSF(SP_Mesh, const std::vector<T> &);

    // Return reference to mesh.
    const OS_Mesh& get_Mesh() const { return *mesh; }

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
// OS_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//
// CCSF explicit constructor

template<class T>
OS_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_) 
    : mesh(mesh_), data(mesh->num_cells()) 
{
    Require (mesh);
    Ensure  (!empty());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
OS_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_, 
		       const std::vector<T> &array)
    : mesh(mesh_), data(array)
{
    Require (mesh)
    Ensure  (data.size() == mesh->num_cells());
    Ensure  (!empty());
}

//===========================================================================//
// class OS_Mesh::CCVF
//
// cell-centered vector fields
//===========================================================================//

template<class T>
class OS_Mesh::CCVF
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
    // SP back to OS_Mesh
    SP_Mesh mesh;
    
    // 2-D field vector, (dimension, num_cells).
    std::vector<std::vector<T> > data;

  public:
    // Inline explicit constructor.
    inline explicit CCVF(SP_Mesh);

    // Additional constructors.
    inline CCVF(SP_Mesh, const std::vector<std::vector<T> > &);

    // Return reference to mesh.
    const OS_Mesh& get_Mesh() const { return *mesh; }

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
// OS_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//
// CCVF explicit constructor

template<class T>
OS_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_)
    : mesh(mesh_), data(mesh->get_Coord().get_dim())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i < mesh->get_Coord().get_dim(); i++)
	data[i].resize(mesh->num_cells());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
OS_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_, 
		       const std::vector<std::vector<T> > &array)
    : mesh(mesh_), data(array)
{
    // check things out
    Ensure (data.size() == mesh->get_Coord().get_dim());
    for (int dim = 0; dim < mesh->get_Coord().get_dim(); dim++)
	Ensure (data[dim].size() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
const T& OS_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim-1][cell-1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
T& OS_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim-1][cell-1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
std::vector<T> OS_Mesh::CCVF<T>::operator()(int cell) const
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
OS_Mesh::CCVF<T>::iterator OS_Mesh::CCVF<T>::begin(int i)
{
    Require(i > 0 && i <= data.size());
    return data[i-1].begin();
}

template<class T>
OS_Mesh::CCVF<T>::const_iterator OS_Mesh::CCVF<T>::begin(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].begin();
}

template<class T>
OS_Mesh::CCVF<T>::iterator OS_Mesh::CCVF<T>::end(int i)
{
    Require(i > 0 && i <= data.size());
    return data[i-1].end();
}

template<class T>
OS_Mesh::CCVF<T>::const_iterator OS_Mesh::CCVF<T>::end(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].end();
}

template<class T>
OS_Mesh::CCVF<T>::size_type OS_Mesh::CCVF<T>::size(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].size();
}

template<class T>
bool OS_Mesh::CCVF<T>::empty(int i) const
{
    Require(i > 0 && i <= data.size());
    return data[i-1].empty();
}

} // end namespace rtt_mc

#endif                          // __mc_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of mc/OS_Mesh.hh
//---------------------------------------------------------------------------//
