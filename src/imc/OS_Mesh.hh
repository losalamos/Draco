//----------------------------------*-C++-*----------------------------------//
// OS_Mesh.hh
// Thomas M. Evans
// Tue Feb  3 16:50:13 1998
//---------------------------------------------------------------------------//
// @> OS_Mesh mesh class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_OS_Mesh_hh__
#define __imctest_OS_Mesh_hh__

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
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include "imctest/Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include <vector>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <cassert>
#include <iostream>
#include <string>

IMCSPACE

// defined namespaces
using std::vector;
using std::fill;
using std::min_element;
using std::max_element;
using std::ostream;
using std::pow;
using std::cout;
using std::endl;
using std::string;

using RNG::Sprng;
    
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
      // data in field
	vector<T> data;

    public:
      // inline explicit constructor
	explicit CCSF(SP<OS_Mesh> mesh_) : mesh(mesh_),
	    data(mesh->num_cells()) { }

      // return reference to mesh
	const OS_Mesh& get_Mesh() const { return *mesh; }

      // subscripting
	T operator()(int cell) const { return data[cell-1]; }
	T& operator()(int cell) { return data[cell-1]; }
    };  

    template<class T>
    class CCVF
    {
    private:

      // SP back to OS_Mesh
	SP<OS_Mesh> mesh;
      // the data array is data(dimension,num_cells) where
      // dimension is 1 (1-D), 2 (2-D), or 3 (3-D)
	vector< vector<T> > data;

    public:
      // inline explicit constructor
	inline explicit CCVF(SP<OS_Mesh>);

      // return reference to mesh
	const OS_Mesh& get_Mesh() const { return *mesh; }

      // subscripting
	inline T operator()(int, int) const;
	inline T& operator()(int, int);
    };  

  // useful typedefs used when working with a mesh
    typedef vector<double> CCSF_a;
    typedef vector< vector<double> > CCVF_a;
    typedef vector<int> CCSF_i;
    typedef vector< vector<int> > CCVF_i;
   
  // temporary typedefs for compiling code until KCC 3.3 is released
    typedef CCSF<double> CCSF_double;
    typedef CCSF<int> CCSF_int;
    typedef CCVF<double> CCVF_double;
    typedef CCVF<int> CCVF_int;

private:
  // base class reference to a derived coord class
    SP<Coord_sys> coord;
  // layout of mesh
    Layout layout;
  // vertices in mesh
    CCVF_a vertex;
  // cell-pairings of cell to its vertices
    CCVF_i cell_pair;
  // area of surfaces on each dimension
    CCVF_a sur;

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
    OS_Mesh(SP<Coord_sys>, Layout &, CCVF_a &, CCVF_i &); 

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

    void print(ostream &, int) const;

  // End_Verbatim 
  // End_Doc 

  // Begin_Doc os_mesh-rint.tex
  // Begin_Verbatim 

  // services required by ALL mesh types used in JAYENNE

  // references to imbedded objects and data required for Parallel_Building
    const Layout& get_Layout() const { return layout; }
    const Coord_sys& get_Coord() const { return *coord; }
    const CCVF_a& get_vertex() const { return vertex; }
    const CCVF_i& get_cell_pair() const { return cell_pair; }

  // required services for transport; 
    int next_cell(int cell, int face) const { return layout(cell, face); }
    int get_cell(const vector<double> &) const;
    double get_db(const vector<double> &, const vector<double> &, int, 
		  int &) const;
    inline vector<double> get_normal(int, int) const;
    inline double volume(int) const;
    vector<int> get_surcells(string) const;
    int get_bndface(string, int) const;
    inline CCVF_a get_vertices(int, int) const;
    inline vector<double> sample_pos(string, int, Sprng &) const;

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

inline ostream& operator<<(ostream &output, const OS_Mesh &object)
{
    for (int cell = 1; cell <= object.num_cells(); cell++)
	object.print(output, cell);
    return output;
}

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// OS_Mesh::CCVF functions
//---------------------------------------------------------------------------//

// CCVF constructor
template<class T>
inline OS_Mesh::CCVF<T>::CCVF(SP<OS_Mesh> mesh_)
    : mesh(mesh_), data(mesh->get_Coord().get_dim())
{
  // initialize data array
    for (int i = 0; i < mesh->get_Coord().get_dim(); i++)
	data[i].resize(mesh->num_cells());
}

//---------------------------------------------------------------------------//

template<class T>
inline T OS_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim-1][cell-1]; 
}

//---------------------------------------------------------------------------//

template<class T>
inline T& OS_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim-1][cell-1];
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

inline vector<double> OS_Mesh::get_normal(int cell, int face) const
{
  // OS_Meshes do not require any functionality from Coord_sys to 
  // calculate the normal, do simple return

  // normal always has 3 components, use Get_sdim()
    vector<double> normal(coord->get_sdim(), 0.0);
	
  // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
  // one coordinate is non-zero    
    normal[(face-1)/2] = pow(-1.0, face);

  // return the normal
    return normal;
}

//---------------------------------------------------------------------------//

inline OS_Mesh::CCVF_a OS_Mesh::get_vertices(int cell, int face) const
{
  // determine the vertices along a cell-face

  // return vertices
    CCVF_a ret_vert(coord->get_dim());

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

  // asserts
    for (int d = 0; d < coord->get_dim(); d++)
	assert (ret_vert[0].size() == ret_vert[d].size());

  // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// sample the position in a cell

inline vector<double> OS_Mesh::sample_pos(string dist, int cell, 
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
    vector<double> r = coord->sample_pos(dist, vmin, vmax, random);

  // return position vector
    return r;
}

CSPACE

#endif                          // __imctest_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Mesh.hh
//---------------------------------------------------------------------------//
