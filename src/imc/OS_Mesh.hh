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
// Date created : 2-3-97
// Purpose      : Orthogonal-Structured Mesh Class
//
// revision history:
// -----------------
// 0) original
// 1) 2-6-98: changed storage of Layout to a full object instead of a smart
//            pointer
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
	    : mesh(mesh_), data(mesh->Layout().getNum_cell(), 0.0),
	      index(mesh->Layout().getNum_cell())
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
	    : mesh(mesh_), data(mesh->Coord().getDim()),
	      index(mesh->Layout().getNum_cell())
	{
	  // initialize data array
	    for (int i = 0; i < mesh->Coord().getDim(); i++)
		data[i].resize(mesh->Layout().getNum_cell());
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
  // index operator for domain-decomposition
    CCF_i index;
  // private copy constructor and assignment operator accessible only by 
  // the OS_Builder
    OS_Mesh(const OS_Mesh &);
    OS_Mesh& operator=(const OS_Mesh &);
public:
  // base class constructor
    OS_Mesh(SP<Coord_sys> coord_, const Layout &layout_, CCVF_a &pos_,
	    CCVF_a &dim_, CCVF_a &sur_, CCF_i &index_)
	: coord(coord_), layout(layout_), pos(pos_), dim(dim_),
	  sur(sur_), index(index_)
    {}
  // member functions
    const Layout& Layout() const { return layout; }
    const Coord_sys& Coord() const { return *coord; }
    double begin(int dim) const { return sur[dim-1].front(); }
    double end(int dim) const { return sur[dim-1].back(); }
    int getCell(vector<double> &) const;
    double getDb(vector<double> &, vector<double> &, int, int &) const;
    void print(int) const;
};

CSPACE

#endif                          // __imctest_OS_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Mesh.hh
//---------------------------------------------------------------------------//
