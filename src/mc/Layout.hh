//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Layout.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 30 15:53:52 1998
 * \brief  Layout class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Layout_hh__
#define __mc_Layout_hh__

#include <iostream>
#include <vector>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class Layout
 *
 * Class that describes the cell-face-cell connectivity for 1-to-1 cell
 * meshes.  1-to-1 describes meshes whose cells each have only 1 cell
 * neighbor across a face.  Layout itself is independent of the geometry of
 * the mesh.  Thus it can be used to describe the connectivity in any
 * structured or unstructered 1-to-1 mesh.
 *
 * This class is used primarily as a component class of valid mesh types.  It
 * is not a general service class.  See rtt_mc::OS_Mesh and
 * rtt_mc::OS_Builder to see how it is used in a mesh implementation.
 */
// revision history:
// -----------------
//  0) original
//  1)  5-20-98 : added overloaded ==,!= operators for desing-by-contract
//                purposes  
// 
//===========================================================================//

class Layout
{
  public:
    // STL Typedefs
    typedef std::vector<std::vector<int> > vf_int;
    
  private:
    // Cell-face-cell array for transport between cells, Layout adjusts the
    // cell and face indices to (cell-1) and (face-1).
    vf_int face_cell;

  public:
    // Default constructor.
    Layout(int num_cells = 0) : face_cell(num_cells) {}

    // Set size member functions.

    // Set size of whole Layout and set size for number of faces in a
    // particular cell.
    void set_size(int num_cells) { face_cell.resize(num_cells); }
    inline void set_size(int, int);

    // Get size member functions.
    int num_cells() const { return face_cell.size(); }
    inline int num_faces(int) const;

    // Diagnostic functions
    void print(std::ostream &, int) const;

    // Overloaded subscripting operators for assignment and retrieval.
    inline int operator()(int, int) const;
    inline int& operator()(int, int);

    // Overloaded operators for equality.
    inline bool operator==(const Layout &) const;
    bool operator!=(const Layout &rhs) const { return !(*this == rhs); }
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// overload operator for stream output

std::ostream& operator<<(std::ostream &, const Layout &);

//---------------------------------------------------------------------------//
// overload equality(==) operator for design-by-contract

bool Layout::operator==(const Layout &rhs) const
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
// set the number of faces for cell cell_index

void Layout::set_size(int cell_index, int num_faces)
{
    face_cell[cell_index-1].resize(num_faces);
} 

//---------------------------------------------------------------------------//
// return the number of faces for cell cell_index

int Layout::num_faces(int cell_index) const
{
    return face_cell[cell_index-1].size();
}

//---------------------------------------------------------------------------//
// subscripting operator for data referencing

int Layout::operator()(int cell_index, int face_index) const
{
    return face_cell[cell_index-1][face_index-1];
}

//---------------------------------------------------------------------------//
// subscripting operator for data assignment

int& Layout::operator()(int cell_index, int face_index)
{
    return face_cell[cell_index-1][face_index-1];
}

} // end namespace rtt_mc

#endif                          // __mc_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Layout.hh
//---------------------------------------------------------------------------//
