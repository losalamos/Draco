//----------------------------------*-C++-*----------------------------------//
// Layout.hh
// B.T. Adams
// 23 Sept 15:53:52 1999
//---------------------------------------------------------------------------//
// @> Layout class header file  - this file should really be moved 
// back to mc. It was modified from the original for AMR meshes to allow a
// cell to have more than one neighbor, but it defaults to a single adjacent
// cell and thus should work just hunky dory with the OS_Mesh type too (with
// a small change to OS_Builder.cc also included to accomodate a ragged right
// array for the adjacent cell vector).
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

#ifndef __mc_Layout_hh__
#define __mc_Layout_hh__

//===========================================================================//
// class Layout - 
//
// Purpose : base class which describes the cell-face-cell interactions for 
//           MC; Layout itself is independent of the geometry of the mesh
//           but Mesh has to build it because Mesh can calculate the neighbor 
//           info from the geometry
//
// revision history:
// -----------------
//  0) original
//  1)  5-20-98 : added overloaded ==,!= operators for desing-by-contract
//                purposes  
// 
//===========================================================================//

#include <iostream>
#include <vector>

namespace rtt_mc 
{

// defined namespaces
using std::ostream;
using std::vector;

class Layout
{
private:
  // cell-face-cell array for transport between cells, Layout adjusts the
  // cell, face, and adjacent cell indices to (cell-1), (face-1), and
  // (adj-1)
    vector<vector<vector<int> > > face_cell;

  // Begin_Doc layout-int.tex
  // Begin_Verbatim 

public:
  // inline default constructor
    Layout(int num_cells = 0) : face_cell(num_cells) {}

  // set size member functions

  // set size of whole Layout, number of faces in a particular cell, and
  // number of adjacent cells.
    void set_size(int num_cells) { face_cell.resize(num_cells); }
    inline void set_size(int num_cells, int num_faces);
    inline void set_adj_size(int cell_index, int face_index, int num_adj = 1);
  // get size member functions
    int num_cells() const { return face_cell.size(); }
    inline int num_faces(int cell) const;
    inline int num_adj(int cell, int face) const;

  // diagnostic functions
    void print(ostream &, int) const;

  // overloaded subscripting operators for assignment and retrieval
    inline int operator()(int cell, int face, int adjcell = 1) const;
    inline int& operator()(int cell, int face, int adjcell = 1);

  // overloaded operators for equality
    inline bool operator==(const Layout &) const;
    bool operator!=(const Layout &rhs) const { return !(*this == rhs); }

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// overload operator for stream output

ostream& operator<<(ostream &, const Layout &);

//---------------------------------------------------------------------------//
// overload equality(==) operator for design-by-contract

inline bool Layout::operator==(const Layout &rhs) const
{
  // if the data is equal, the Layouts are equal
    if (face_cell == rhs.face_cell)
	return true;
    
  // if we haven't returned then the Layouts aren't equal
    return false;
}

//---------------------------------------------------------------------------//
// inline functions for Layout
//---------------------------------------------------------------------------//
// set the number of faces for cell cell_index, and the 

inline void Layout::set_size(int cell_index, int num_faces)
{
    face_cell[cell_index-1].resize(num_faces);
} 

//---------------------------------------------------------------------------//
// set the number of adjcacent cells for face face_index of cell cell_index

inline void Layout::set_adj_size(int cell_index, int face_index, int num_adj)
{
    face_cell[cell_index-1][face_index-1].resize(num_adj);
} 

//---------------------------------------------------------------------------//
// return the number of faces for cell cell_index

inline int Layout::num_faces(int cell_index) const
{
    return face_cell[cell_index-1].size();
}

//---------------------------------------------------------------------------//
// return the number of adjacent cells for face face_index of cell cell_index

inline int Layout::num_adj(int cell_index, int face_index) const
{
    return face_cell[cell_index-1][face_index-1].size();
}

//---------------------------------------------------------------------------//
// subscripting operator for data referencing

inline int Layout::operator()(int cell_index, int face_index, 
			      int adj_index) const
{
    return face_cell[cell_index-1][face_index-1][adj_index-1];
}

//---------------------------------------------------------------------------//
// subscripting operator for data assignment

inline int& Layout::operator()(int cell_index, int face_index, 
			       int adj_index)
{
    return face_cell[cell_index-1][face_index-1][adj_index-1];
}

} // end namespace rtt_mc

#endif                          // __mc_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Layout.hh
//---------------------------------------------------------------------------//
