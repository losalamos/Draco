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
// Purpose : base class which describes the cell-face-cell inter-
//           actions for IMC; Layout itself is independent of
//           the geometry of the mesh but Mesh has to build it
//           because Mesh can calculate the neighbor info from
//           the geometry
//
// revision history:
// -----------------
// 0) original
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
public:
  // inline default constructor, can give the total number of
  // cells as an argument; no copy constructor or assignment
  // operator needed because we are using vectors
    Layout(int num_cells = 0) : face_cell(num_cells) {}

  // inline function which sets the size of the face_cell
  // array to the number of cells in the problem 
    void Set_size(int num_cells) { face_cell.resize(num_cells); }

  // inline function to size the face vector for cell_index
    void Set_size(int cell_index, int num_faces)
    {
	face_cell[cell_index-1].resize(num_faces);
    } 

  // inline get functions
    int Num_cells() const { return face_cell.size(); }
    int Num_faces(int cell_index) const
    {
	return face_cell[cell_index-1].size();
    }
    void Print(int) const;

  // overloaded operator for subscripting, not assignment,
  // for const objects
    int operator()(int cell_index, int face_index) const
    {
	return face_cell[cell_index-1][face_index-1];
    }

  // overloaded operator for subscriping and assignment
    int & operator()(int cell_index, int face_index)
    {
	return face_cell[cell_index-1][face_index-1];
    }
};

// overload operator declarations for Layout
ostream & operator<<(ostream &, const Layout &);

CSPACE

#endif                          // __imctest_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Layout.hh
//---------------------------------------------------------------------------//
