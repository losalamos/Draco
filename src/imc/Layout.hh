//----------------------------------*-C++-*----------------------------------//
// Layout.hh
// Thomas M. Evans
// Fri Jan 30 15:53:52 1998
//---------------------------------------------------------------------------//
// @> Layout class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Layout_hh__
#define __imctest_Layout_hh__

//===========================================================================//
// class Layout - 
//
// Purpose : base class which describes the cell-face-cell interactions for 
//           IMC; Layout itself is independent of the geometry of the mesh
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

#include "imctest/Names.hh"
#include <iostream>
#include <vector>

IMCSPACE

// defined namespaces
using std::ostream;
using std::vector;

class Layout
{
private:
  // cell-face-cell array for transport between cells, Layout adjusts the
  // cell and face indices to (cell-1) and (face-1)
    vector< vector<int> > face_cell;

  // Begin_Doc layout-int.tex
  // Begin_Verbatim 

public:
  // inline default constructor
    Layout(int num_cells = 0) : face_cell(num_cells) {}

  // set size member functions

  // set size of whole Layout and set size for number of faces in a
  // particular cell
    void set_size(int num_cells) { face_cell.resize(num_cells); }
    inline void set_size(int, int);

  // get size member functions
    int num_cells() const { return face_cell.size(); }
    inline int num_faces(int) const;

  // diagnostic functions
    void print(ostream &, int) const;

  // overloaded subscripting operators for assignment and retrieval
    inline int operator()(int, int) const;
    inline int& operator()(int, int);

  // overloaded operators for equality
    inline bool operator==(const Layout &) const;
    bool operator!=(const Layout &rhs) const { return !(*this == rhs); }

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// overload operator declarations for Layout

ostream& operator<<(ostream &, const Layout &);

//---------------------------------------------------------------------------//
// overload equality(==) operator for design-by-contract

inline bool Layout::operator==(const Layout &rhs) const
{
  // if the data is equal, the Layouts are equal
    if ( face_cell == rhs.face_cell)
	return true;
    
  // if we haven't returned then the Layouts aren't equal
    return false;
}

//---------------------------------------------------------------------------//
// inline functions for Layout
//---------------------------------------------------------------------------//
// set the number of faces for cell cell_index

inline void Layout::set_size(int cell_index, int num_faces)
{
    face_cell[cell_index-1].resize(num_faces);
} 

//---------------------------------------------------------------------------//
// return the number of faces for cell cell_index

inline int Layout::num_faces(int cell_index) const
{
    return face_cell[cell_index-1].size();
}

//---------------------------------------------------------------------------//
// subscripting operator for data referencing

inline int Layout::operator()(int cell_index, int face_index) const
{
    return face_cell[cell_index-1][face_index-1];
}

//---------------------------------------------------------------------------//
// subscripting operator for data assignment

inline int& Layout::operator()(int cell_index, int face_index)
{
    return face_cell[cell_index-1][face_index-1];
}

CSPACE

#endif                          // __imctest_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Layout.hh
//---------------------------------------------------------------------------//
