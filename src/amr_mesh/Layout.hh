//----------------------------------*-C++-*----------------------------------//
// Layout.hh
// B.T. Adams
// 23 Sept 15:53:52 1999
/*! 
 * \file   amr_mesh/Layout.hh
 * \author B.T. Adams
 * \date   Wed Sep 15 10:33:26 1999
 * \brief  Header file for Layout class library.
 */
//---------------------------------------------------------------------------//
// @> Layout class header file  - this file should really be moved 
// back to mc. It was modified from the original for AMR meshes to allow a
// cell to have more than one neighbor, but it defaults to a single adjacent
// cell and thus should work just hunky dory with the OS_Mesh type too (with
// a small change to OS_Builder.cc also included to accomodate a ragged right
// array for the adjacent cell vector).
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

#ifndef __amr_Layout_hh__
#define __amr_Layout_hh__

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
//  1)  5-20-98 : added overloaded ==,!= operators for design-by-contract
//                purposes  
// 
//===========================================================================//

#include <iostream>
#include <vector>

namespace rtt_amr 
{

// defined namespaces
using std::ostream;
using std::vector;

/*!
 * \brief  Layout is a base class that describes the relation of a cell's 
 *         faces to those of the neighboring cells (i.e., the "connectivity").
 *         The Layout class is independent of the geometry of the mesh, but a 
 *         Mesh-type class requires the Layout to calculate neighbor 
 *         information from the geometry data.
 *
 *\sa The Layout class is used by both the OS_Mesh and CAR_CU_Mesh classes. 
 *    An \ref amr_overview is provided to describe the basic functionality 
 *    of that particular mesh class (which is also very similar to the 
 *    functionality of the OS_Mesh class from which it was derived).
 */     
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
/*!
 * \brief Constructs a Layout class object.
 * \param  num_cells Number of cells in the mesh (defaults to zero). 
 */
    Layout(int num_cells = 0) : face_cell(num_cells) {}

  // set size member functions

  // set size of whole Layout, number of faces in a particular cell, and
  // number of adjacent cells.
/*!
 * \brief Resizes the Layout class object to accomodate the specified number 
 *        of cells in the mesh
 * \param  num_cells Number of cells in the mesh. 
 */
    void set_size(int num_cells) { face_cell.resize(num_cells); }
/*!
 * \brief Resizes the specified cell of the Layout class object to accomodate
 *        the specified number of faces per cell.
 * \param  cell_index Cell number. 
 * \param  num_faces Number of faces per cell. 
 */
    inline void set_size(int cell_index, int num_faces);
/*!
 * \brief Resizes the specified cell face of the Layout class object to 
 *        accomodate the specified number of adjacent cells.
 * \param  cell_index Cell number. 
 * \param  face_index Cell face number. 
 * \param  num_adj Number of adjacent cells.
 */
    inline void set_adj_size(int cell_index, int face_index, int num_adj = 1);
  // get size member functions
/*!
 * \brief Returns the number of cells in the mesh.
 * \return Number of mesh cells. 
 */
    int num_cells() const { return face_cell.size(); }
/*!
 * \brief Returns the number of faces for the specified cell.
 * \param  cell Cell number. 
 * \return Number of cell faces. 
 */
    inline int num_faces(int cell) const;
/*!
 * \brief Returns the number of cells adjacent to the specified cell face.
 * \param  cell Cell number. 
 * \param  face Cell face number. 
 * \return Number of cells adjacent to the cell face. 
 */
    inline int num_adj(int cell, int face) const;

  // diagnostic functions
/*!
 * \brief Diagnostic member function used to print out the cells adjacent to
 *        each face of the specified cell.
 * \param output Stream-output class object.
 * \param  cell_index Cell number.
 */
    void print(ostream & output, int cell_index) const;

  // overloaded subscripting operators for assignment and retrieval
/*!
 * \brief Overloaded operator to return the cell number of the specified 
 *        adjacent cell index for the specified cell face.
 * \param cell Cell number.
 * \param face Cell face number. 
 * \param adjcell Adjacent cell index (defaults to 1).
 * \return Adjacent cell number.
 */
    inline int operator()(int cell, int face, int adjcell = 1) const;
/*!
 * \brief Overloaded operator to assign a cell number to the specified 
 *        adjacent cell index for the specified cell face.
 * \param cell Cell number.
 * \param face Cell face number. 
 * \param adjcell Adjacent cell index (defaults to 1).
 * \return Adjacent cell number.
 */
    inline int & operator()(int cell, int face, int adjcell = 1);

  // overloaded operators for equality
/*!
 * \brief Compares a layout for equivalence relative to the current layout.
 * \param rhs Second layout
 * \return Status of second layout = current layout. 
 */
    inline bool operator==(const Layout & rhs) const;
/*!
 * \brief Compares a layout for non-equivalence relative to the current layout.
 * \param rhs Second layout
 * \return Status of second layout != current layout. 
 */
    bool operator!=(const Layout & rhs) const { return !(*this == rhs); }

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// overload operator for stream output

/*!
 * \brief Overloaded stream-insertion operator for layout output.
 * \param output Stream-output class object.
 * \param object Layout class object.
 * \return Reference to output.
 */
ostream & operator<<(ostream & output, const Layout & object)
{
    int num_cells = object.num_cells();
    for (int cell_index = 1; cell_index <= num_cells; cell_index++)
        object.print(output, cell_index);
    return output;
}

//---------------------------------------------------------------------------//
// overload equality(==) operator for design-by-contract

inline bool Layout::operator==(const Layout & rhs) const
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
// set the number of adjacent cells for face face_index of cell cell_index

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

} // end namespace rtt_amr

#endif                          // __amr_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Layout.hh
//---------------------------------------------------------------------------//
