//----------------------------------*-C++-*----------------------------------//
// OS_Mesh.hh
// Thomas M. Evans
// Tue Feb  3 16:50:13 1998
//---------------------------------------------------------------------------//
// @> OS_Mesh mesh class header file
//---------------------------------------------------------------------------//

#ifndef __imc_OS_Mesh_hh__
#define __imc_OS_Mesh_hh__

//===========================================================================//
// class OS_Mesh - 
//
// Purpose : Orthogonal-Structured Mesh Class
//
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
// 
//===========================================================================//

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

namespace rtt_imc 
{

// stl namespaces
using std::vector;
using std::fill;
using std::min_element;
using std::max_element;
using std::ostream;
using std::pow;
using std::endl;
using std::string;

// draco namespaces
using rtt_rng::Sprng;
using dsxx::SP;
    
class OS_Mesh
{
public:
  // The cell centered fields are used only outside of the MESH-type class

  // class definitions of the cell-centered fields: neither of these classes
  // require copy constructors or assignment operators as the SP<> and 
  // vector<> classes can do assignment
    template<class T>
    class CCSF
    {
    private:
      // SP back to OS_Mesh 
	SP<OS_Mesh> mesh;
      // data in field, (num_cells)
	vector<T> data;

    public:
      // inline explicit constructor
	inline explicit CCSF(SP<OS_Mesh>);

      // additional constructors
	inline CCSF(SP<OS_Mesh>, const vector<T> &);

      // return reference to mesh
	const OS_Mesh& get_Mesh() const { return *mesh; }

      // subscripting
	const T& operator()(int cell) const { return data[cell-1]; }
	T& operator()(int cell) { return data[cell-1]; }
    };  

    template<class T>
    class CCVF
    {
    private:

      // SP back to OS_Mesh
	SP<OS_Mesh> mesh;
      // 2-D field vector, (dimension, num_cells)
	vector< vector<T> > data;

    public:
      // inline explicit constructor
	inline explicit CCVF(SP<OS_Mesh>);

      // additional constructors
	inline CCVF(SP<OS_Mesh>, const vector<vector<T> > &);

      // return reference to mesh
	const OS_Mesh& get_Mesh() const { return *mesh; }

      // subscripting
	inline const T& operator()(int, int) const;
	inline T& operator()(int, int);

      // getting a CC vector
	inline vector<T> operator()(int) const;
    };  

  // useful typedefs used when working with a mesh
    typedef vector<double> CCSF_a;
    typedef vector< vector<double> > CCVF_d;
    typedef vector<int> CCSF_i;
    typedef vector< vector<int> > CCVF_i;
   
  // temporary typedefs for compiling code until KCC 3.3+ is released
    typedef CCSF<double> CCSF_double;
    typedef CCSF<int> CCSF_int;
    typedef CCVF<double> CCVF_double;
    typedef CCVF<int> CCVF_int;
    typedef CCSF<string> CCSF_string;

private:
  // base class reference to a derived coord class
    SP<Coord_sys> coord;
  // layout of mesh
    Layout layout;
  // vertices in mesh
    CCVF_d vertex;
  // cell-pairings of cell to its vertices
    CCVF_i cell_pair;
  // area of surfaces on each dimension
    CCVF_d sur;
  // indicator whether this is a submesh
    bool submesh;

  // private functions

  // calculate a surface array from the vertices of the mesh
    void calc_surface();

  // private copy and assignment operators; can't copy or assign a mesh
    OS_Mesh(const OS_Mesh &);
    OS_Mesh& operator=(const OS_Mesh &);

  // Begin_Doc os_mesh-int.tex
  // Begin_Verbatim 

public:
  // generalized constructor for all mesh types
    OS_Mesh(SP<Coord_sys>, Layout &, CCVF_d &, CCVF_i &, bool = false); 

  // member functions used by the OS_Mesh-dependent classes

  // mesh dimensionality functions

    inline double begin(int) const;
    inline double end(int) const;
    int num_cells() const { return layout.num_cells(); }

  // cell dimensionality functions

  // find minimum and maximum dimension of cell
    inline double min(int, int) const;
    inline double max(int, int) const;

  // find centerpoint of cell and width of cell
    inline double pos(int, int) const;
    double dim(int d, int cell) const { return max(d, cell) - min(d, cell); } 

  // diagnostic functions
    void print(ostream &) const;
    void print(ostream &, int) const;

  // End_Verbatim 
  // End_Doc 

  // Begin_Doc os_mesh-rint.tex
  // Begin_Verbatim 

  // services required by ALL mesh types used in JAYENNE

  // references to imbedded objects and data required for Parallel_Building
    const Layout& get_Layout() const { return layout; }
    const Coord_sys& get_Coord() const { return *coord; }
    SP<Coord_sys> get_SPCoord() const { return coord; }
    const CCVF_d& get_vertex() const { return vertex; }
    const CCVF_i& get_cell_pair() const { return cell_pair; }

  // required services for transport; 
    int next_cell(int cell, int face) const { return layout(cell, face); }
    int get_cell(const vector<double> &) const;
    double get_db(const vector<double> &, const vector<double> &, int, 
		  int &) const;
    inline vector<double> get_normal(int, int) const;
    inline vector<double> get_normal_in(int, int) const;
    inline double volume(int) const;
    inline double face_area(int, int) const;
    vector<int> get_surcells(string) const;
    int get_bndface(string, int) const;
    inline CCVF_d get_vertices(int, int) const;
    inline CCVF_d get_vertices(int) const;
    inline vector<double> sample_pos(int, Sprng &) const;
    inline vector<double> sample_pos(int, Sprng &, vector<double>, 
				     double) const; 
    inline vector<double> sample_pos_on_face(int, int, Sprng &)	const; 

  // overloaded operators
    bool operator==(const OS_Mesh &) const;
    bool operator!=(const OS_Mesh &rhs) const { return !(*this == rhs); }

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

inline ostream& operator<<(ostream &output, const OS_Mesh &object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// OS_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//
// CCSF explicit constructor

template<class T>
inline OS_Mesh::CCSF<T>::CCSF(SP<OS_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->num_cells()) 
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline OS_Mesh::CCSF<T>::CCSF(SP<OS_Mesh> mesh_, const vector<T> &array)
    : mesh(mesh_), data(array)
{
  // make sure things are kosher
    Ensure (data.size() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// OS_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//
// CCVF explicit constructor

template<class T>
inline OS_Mesh::CCVF<T>::CCVF(SP<OS_Mesh> mesh_)
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
inline OS_Mesh::CCVF<T>::CCVF(SP<OS_Mesh> mesh_, 
			      const vector<vector<T> > &array)
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
inline const T& OS_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim-1][cell-1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T& OS_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim-1][cell-1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
inline vector<T> OS_Mesh::CCVF<T>::operator()(int cell) const
{
  // declare return vector
    vector<T> x;
    
  // loop through dimensions and make return vector for this cell
    for (int i = 0; i < data.size(); i++)
	x.push_back(data[i][cell-1]);

  // return
    Ensure (x.size() == data.size());
    return x;
}

//---------------------------------------------------------------------------//
// OS_Mesh inline functions
//---------------------------------------------------------------------------//

inline double OS_Mesh::begin(int d) const 
{
  // find the minimum surface for d over the whole mesh
    return *min_element(vertex[d-1].begin(), vertex[d-1].end()); 
}

//---------------------------------------------------------------------------//

inline double OS_Mesh::end(int d) const 
{
  // find the maximum surface for d over the whole mesh
    return *max_element(vertex[d-1].begin(), vertex[d-1].end()); 
}

//---------------------------------------------------------------------------//

inline double OS_Mesh::pos(int d, int cell) const
{
  // find center position of cell

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

inline double OS_Mesh::min(int d, int cell) const 
{	
  // find minimum dimension along d of cell

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

inline double OS_Mesh::max(int d, int cell) const
{
  // find maximum dimension of cell

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

inline double OS_Mesh::volume(int cell) const 
{
  // calculate volume of cell

  // loop through dimensions and get volume
    double volume = 1.0;
    for (int d = 1; d <= coord->get_dim(); d++)
	volume *= dim(d, cell);

  // return volume
    return volume;
}

//---------------------------------------------------------------------------//

inline double OS_Mesh::face_area(int cell, int face) const 
{
  // calculate area of face on cell

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

inline vector<double> OS_Mesh::get_normal(int cell, int face) const
{
  // OS_Meshes do not require any functionality from Coord_sys to 
  // calculate the outward normal, do simple return

  // normal always has 3 components, use Get_sdim()
    vector<double> normal(coord->get_sdim(), 0.0);
	
  // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
  // one coordinate is non-zero    
    normal[(face-1)/2] = pow(-1.0, face);

  // return the normal
    return normal;
}

//---------------------------------------------------------------------------//

inline vector<double> OS_Mesh::get_normal_in(int cell, int face) const
{
  // OS_Meshes do not require any functionality from Coord_sys to 
  // calculate the inward normal, do simple return

  // normal always has 3 components, use Get_sdim()
    vector<double> normal(coord->get_sdim(), 0.0);
	
  // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
  // one coordinate is non-zero    
    normal[(face-1)/2] = pow(-1.0, face-1);

  // return the normal
    return normal;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell face

inline OS_Mesh::CCVF_d OS_Mesh::get_vertices(int cell, int face) const
{
  // determine the vertices along a cell-face

  // return vertices
    CCVF_d ret_vert(coord->get_dim());

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

inline OS_Mesh::CCVF_d OS_Mesh::get_vertices(int cell) const
{
  // determine the vertices bounding a cell
    
  // return vertices
    CCVF_d ret_vert(coord->get_dim());

  // loop over cell vertices and build the cell vertices
    for (int i = 0; i < cell_pair[cell-1].size(); i++)
	for (int d = 0; d < coord->get_dim(); d++)
	    ret_vert[d].push_back(vertex[d][cell_pair[cell-1][i]-1]);

  // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// sample the position uniformly in a cell

inline vector<double> OS_Mesh::sample_pos(int cell, Sprng &random) const
{
  // assign minimums and maximums for cell dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

  // use coord_sys to sample the location
    vector<double> r = coord->sample_pos(vmin, vmax, random);

  // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample the position in a cell given a tilt (or other slope-like function)

inline vector<double> OS_Mesh::sample_pos(int cell, Sprng &random, 
					  vector<double> slope, 
					  double center_pt) const
{
  // assign minimums and maximums for cells dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

  // use coord_sys to sample the location
    vector<double> r = coord->sample_pos(vmin, vmax, random, slope,
					 center_pt);

  // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample a position on a face

inline vector<double> OS_Mesh::sample_pos_on_face(int cell, int face, 
						  Sprng &random) const
{
  // assign minimums and maximums for cell dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

  // use coord_sys to sample the location
    vector<double> r = coord->sample_pos_on_face(vmin, vmax, face, random);

  // return position vector
    return r;
}

} // end namespace rtt_imc

#endif                          // __imc_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imc/OS_Mesh.hh
//---------------------------------------------------------------------------//
