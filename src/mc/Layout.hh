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

#include "ds++/SP.hh"

#include <iostream>
#include <vector>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class Layout

 * \brief Provides cell connectivity information for 1-1 meshes.

 * Class that describes the cell-face-cell connectivity for 1-to-1 cell
 * meshes.  1-to-1 describes meshes whose cells each have only 1 cell
 * neighbor across a face.  Layout itself is independent of the geometry of
 * the mesh.  Thus it can be used to describe the connectivity in any
 * structured or unstructered 1-to-1 mesh.

 * This class is used primarily as a component class of valid mesh types.  It
 * is not a general service class.  See rtt_mc::OS_Mesh and
 * rtt_mc::OS_Builder to see how it is used in a mesh implementation.

 * There is no "defined" face ordering.  One can build the Layout however one
 * desires.  However, the face numbering must be consistent throughout the
 * Layout to get proper performance.

 * Data is built into Layout through the use of overloaded operator().

 */
// revision history:
// -----------------
//  0) original
//  1)  5-20-98 : added overloaded ==,!= operators for design-by-contract
//                purposes  
//  2) 04-18-01 : added packing capability
// 
//===========================================================================//

class Layout
{
  public:
    // Forward declaration of pack class.
    struct Pack;

  public:
    // STL Typedefs
    typedef std::vector<int>               sf_int;
    typedef std::vector<std::vector<int> > vf_int;
    typedef rtt_dsxx::SP<Layout::Pack>     SP_Pack;
    typedef rtt_dsxx::SP<Layout>           SP_Layout;
    
  private:
    // Cell-face-cell array for transport between cells, Layout adjusts the
    // cell and face indices to (cell-1) and (face-1).
    vf_int face_cell;

  public:
    // Default constructor.
    explicit Layout(int num_cells = 0);

    // Constructor when the number of cells per face are known.
    Layout(int, int);

    // >>> ACCESSORS

    //! Get the number of cells in the layout.
    int num_cells() const { return face_cell.size(); }

    // Get the number of faces for a given cell.
    inline int num_faces(int) const;

    // >>> MODIFIERS

    //! Set size of layout to number of cells.
    void set_size(int num_cells) { face_cell.resize(num_cells); }

    // Set the number of faces for a given cell.
    inline void set_size(int, int);

    // >>> OVERLOADED OPERATORS

    // Overloaded subscripting operators for assignment and retrieval.
    inline int  operator()(int, int) const;
    inline int& operator()(int, int);

    // Pack function.
    SP_Pack pack(const sf_int &) const;

    // Overloaded operators for equality.
    bool operator==(const Layout &) const;
    bool operator!=(const Layout &rhs) const { return !(*this == rhs); }

    // Diagnostic functions.
    void print(std::ostream &, int) const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
//! Overloaded operator for stream output.

std::ostream& operator<<(std::ostream &, const Layout &);

//---------------------------------------------------------------------------//
// LAYOUT INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*! 
 * \brief Set the number of faces in a cell.

 * \param cell_index cell that faces are being set
 * \param num_faces number of faces to set for cell_index

 */
void Layout::set_size(int cell_index, int num_faces)
{
    face_cell[cell_index-1].resize(num_faces);
} 

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of faces for cell.
 */
int Layout::num_faces(int cell_index) const
{
    return face_cell[cell_index-1].size();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the cell index across a face.

 * Returns the cell index of the cell adjacent to cell_index and across the
 * face given by face_index.

 * This operator only allows access.

 */
int Layout::operator()(int cell_index, int face_index) const
{
    return face_cell[cell_index-1][face_index-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the cell index across a face.

 * Returns the cell index of the cell adjacent to cell_index and across the
 * face given by face_index.

 * This operator may be used as an L-value to change the cell across the
 * face.

 */
int& Layout::operator()(int cell_index, int face_index)
{
    return face_cell[cell_index-1][face_index-1];
}

//===========================================================================//
/*!  
 * \struct Layout::Pack 
 
 * \brief Nested class for packing layouts into raw data for writing
 * or communication.
 
 */
//===========================================================================//

struct Layout::Pack
{
  private:
    // Data contained in the Layout.
    int *data;
    int  size;

    // Disallow assignment.
    const Pack& operator=(const Pack &);

  public:
    // Constructor.
    Pack(int, int *);
    
    // Copy constructor.
    Pack(const Pack &);

    // Destructor.
    ~Pack();

    // >>> Accessors.
    
    //! Get pointer to beginning of integer data stream.
    const int* begin() const { return &data[0]; }
    
    //! Get pointer to end of integer data stream.
    const int* end() const { return &data[size]; }

    //! Get size of integer data stream.
    int get_size() const { return size; }

    //! Get number of packed cells in the layout.
    int get_num_packed_cells() const { return data[0]; }

    // Unpack function.
    SP_Layout unpack() const;
};

} // end namespace rtt_mc

#endif                          // __mc_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Layout.hh
//---------------------------------------------------------------------------//
