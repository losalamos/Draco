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
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include "imctest/Layout.hh"
#include "ds++/SP.hh"
#include <vector>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <cassert>
#include <iostream>

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
    
class OS_Mesh
{
public:
  // The cell centered fields are used only outside of the MESH-type class
  // typedefs to container classes used in the CC classes
    typedef vector<double> CCSF_a;
    typedef vector< vector<double> > CCVF_a;
    typedef vector<int> CCSF_i;
    typedef vector< vector<int> > CCVF_i;

  // class definitions of the cell-centered fields: neither of these classes
  // require copy constructors or assignment operators as the SP<> and 
  // vector<> classes can do assignment
    template<class T>
    class CCSF
    {
    private:
      // SP back to OS_Mesh for getting OS_Mesh stuff like
      // number of cells, etc.
	SP<OS_Mesh> mesh;
      // data in field
	vector<T> data;
    public:
      // inline explicit constructor, must give a OS_Mesh,
      // no copy or operator+ needed
	explicit CCSF(SP<OS_Mesh> mesh_)
	    : mesh(mesh_), data(mesh->Num_cells())
	{ }

      // return reference to mesh
	const OS_Mesh& Mesh() const { return *mesh; }

      // subscripting
	T operator()(int cell) const { return data[cell-1]; }
	T& operator()(int cell_index) { return data[cell_index-1]; }
    };  

    template<class T>
    class CCVF
    {
    private:
      // SP back to OS_Mesh, like CCSF, a CCVF cannot
      // be defined without a mesh
	SP<OS_Mesh> mesh;
      // the data array is data(dimension,num_cells) where
      // dimension is 1 (1-D), 2 (2-D), or 3 (3-D)
	vector< vector<T> > data;
    public:
      // inline explicit constructor, must give a OS_Mesh,
      // no copy or assignment operator needed
	explicit CCVF(SP<OS_Mesh> mesh_)
	    : mesh(mesh_), data(mesh->Coord().Get_dim())
	{
	  // initialize data array
	    for (int i = 0; i < mesh->Coord().Get_dim(); i++)
		data[i].resize(mesh->Num_cells());
	}

      // return reference to mesh
	const OS_Mesh& Mesh() const { return *mesh; }

      // subscripting
	T operator()(int dim, int cell) const 
	{
	    return data[dim-1][cell-1]; 
	}
	T& operator()(int dim, int cell)
	{
	    return data[dim-1][cell-1];
	}
    };  

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
    void Calc_surface();

  // private copy and assignment operators; can't copy or assign a mesh
    OS_Mesh(const OS_Mesh &);
    OS_Mesh& operator=(const OS_Mesh &);
public:
  // default constructor
    OS_Mesh(SP<Coord_sys>, Layout &, CCVF_a &, CCVF_i &); 

  // member functions

  // references to imbedded objects

  // return references to the Coord_sys and Layout
    const Layout& Get_Layout() const { return layout; }
    const Coord_sys& Coord() const { return *coord; }

  // mesh dimensionality functions
    double Begin(int d) const 
    { 
	return *min_element(vertex[d-1].begin(), vertex[d-1].end()); 
    }

    double End(int d) const 
    {
	return *max_element(vertex[d-1].begin(), vertex[d-1].end()); 
    }

    int Num_cells() const { return layout.Num_cells(); }

  // cell dimensionality functions

  // find centerpoint of cell
    double Pos(int d, int cell) const
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

  // find minimum dimension of cell
    double Min(int d, int cell) const 
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

  // find maximum dimension of cell
    double Max(int d, int cell) const
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

  // find the width along dimension of cell
    double Dim(int d, int cell) const { return Max(d, cell) - Min(d, cell); }

  // calculate volume of cell
    double Volume(int cell) const 
    {
      // loop through dimensions and get volume
	double volume = 1.0;
	for (int d = 1; d <= coord.Get_dim(); d++)
	    volume *= Dim(d, cell);

      // return volume
	return volume;
    }

  // diagnostic functions
    void Print(ostream &, int) const;

  // required services
    int Next_cell(int cell, int face) const { return layout(cell, face); }
    int Get_cell(const vector<double> &) const;
    double Get_db(const vector<double> &, const vector<double> &, int, 
		  int &) const;
    vector<double> Get_normal(int cell, int face) const
    {
      // OS_Meshes do not require any functionality from Coord_sys to 
      // calculate the normal, do simple return

      // normal always has 3 components, use Get_sdim()
	vector<double> normal(coord->Get_sdim(), 0.0);
	
      // calculate normal based on face, (-x, +x, -y, +y, -z, +z), only
      // one coordinate is non-zero
	normal[(face-1)/2] = pow(-1.0, face);

      // return the normal
	return normal;
    }
};

CSPACE

#endif                          // __imctest_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Mesh.hh
//---------------------------------------------------------------------------//
