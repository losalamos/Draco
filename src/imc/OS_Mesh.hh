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
// 
//===========================================================================//

// we must include the definitions to Layout and Coord_sys because
// the OS_Mesh class allows access to member functions of the class,
// ie. a class which HoldsA/UsesA OS_Mesh accesses Layout and Coord_sys
// through mesh; also, because CCSF and CCVF are all inline functions,
// it needs to know the member functions of Layout and Coord_sys
#include "Names.hh"
#include "Coord_sys.hh"
#include "Layout.hh"
#include "SP.hh"
#include <vector>
#include <algorithm>
#include <cassert>

IMCSPACE

// defined namespaces
using std::vector;
using std::fill;
    
class OS_Mesh
{
public:
  // The cell centered fields are used only outside of the MESH-type class
  // typedefs to container classes used in the CC classes
    typedef vector<double> CCSF_a;
    typedef vector< vector<double> > CCVF_a;
    typedef vector<int> CCF_i;
  // OS Builder is a friend so that it may use the private copy constructor
  // and assignment operators
    friend class OS_Builder;

  // class definitions of the cell-centered fields
    class CCSF
    {
    private:
      // SP back to OS_Mesh for getting OS_Mesh stuff like
      // number of cells, etc.
	SP<OS_Mesh> mesh;
      // data in field
	CCSF_a data;
	CCF_i index;
    public:
      // inline explicit constructor, must give a OS_Mesh,
      // no copy or operator+ needed
	explicit CCSF(SP<OS_Mesh> mesh_)
	    : mesh(mesh_), data(mesh->Num_cells(), 0.0),
	      index(mesh->Num_cells())
	{
	  // initialize index array to zero, remember, cannot initialize
	  // a vector<int> with two arguments of the same type because
	  // the iterator constructor is called instead
	    fill(index.begin(), index.end(), 0);
	}
      // allow OS_Mesh info to be accessed from CCSF
	const OS_Mesh& Mesh() const { return *mesh; }
      // overloaded operator for subscripting, not assignment
	double operator()(int cell_index) const
	{
	    return data[cell_index-1];
	}
      // overloaded operator for subscripting and assignment
	double& operator()(int cell_index)
	{
	    return data[cell_index-1];
	}
    };  

    class CCVF
    {
    private:
      // SP back to OS_Mesh, like CCSF, a CCVF cannot
      // be defined without a mesh
	SP<OS_Mesh> mesh;
      // the data array is data(dimension,num_cells) where
      // dimension is 1 (1-D), 2 (2-D), or 3 (3-D)
	CCVF_a data;
	CCF_i index;
    public:
      // inline explicit constructor, must give a OS_Mesh,
      // no copy or assignment operator needed
	explicit CCVF(SP<OS_Mesh> mesh_)
	    : mesh(mesh_), data(mesh->Coord().Get_dim()),
	      index(mesh->Num_cells())
	{
	  // initialize data array
	    for (int i = 0; i < mesh->Coord().Get_dim(); i++)
		data[i].resize(mesh->Num_cells());
	  // initialize index array to zero
	    fill(index.begin(), index.end(), 0);
	}
	const OS_Mesh& Mesh() const { return *mesh; }
      // overloaded operator for subscripting, not assignment
	double operator()(int dim, int cell_index) const
	{
	    return data[dim-1][cell_index-1];
	}
      // overloaded operator for subscripting and assignment
	double& operator()(int dim, int cell_index)
	{
	    return data[dim-1][cell_index-1];
	}
    };

private:
  // base class reference to a derived coord class
    SP<Coord_sys> coord;
  // layout of mesh
    Layout layout;
  // origin of each cell
    CCVF_a pos;
  // dr dimensions of each cell orthogonal through origin
    CCVF_a dim;
  // surfaces along each dimension axis
    CCVF_a sur;
  // private copy constructor and assignment operator accessible only by 
  // the OS_Builder
    OS_Mesh(const OS_Mesh &);
    OS_Mesh& operator=(const OS_Mesh &);
public:
  // base class constructor
    OS_Mesh(SP<Coord_sys> coord_, const Layout &layout_, CCVF_a &pos_,
	    CCVF_a &dim_, CCVF_a &sur_)
	: coord(coord_), layout(layout_), pos(pos_), dim(dim_), sur(sur_)
    {
	int num_cells = Num_cells();
	assert (num_cells == pos.size());
	assert (num_cells == dim.size());
	assert (num_cells == (sur.size() - 1));
    }

  // member functions

  // helper functions
    const Layout& Layout() const { return layout; }
    const Coord_sys& Coord() const { return *coord; }
    double Begin(int d) const { return sur[d-1].front(); }
    double End(int d) const { return sur[d-1].back(); }
    double Pos(int d, int cell) const { return pos[d-1][cell-1]; }
    double Dim(int d, int cell) const { return dim[d-1][cell-1]; }
    int Num_cells() const { return layout.Num_cells(); }

  // diagnostic functions
    void Print(int) const;

  // required services
    int Next_cell(int cell, int face) const { return layout(cell, face); }
    int Get_cell(vector<double> &) const;
    double Get_db(vector<double> &, vector<double> &, int, int &) const;

};

CSPACE

#endif                          // __imctest_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Mesh.hh
//---------------------------------------------------------------------------//
