//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/AMR_Layout.hh
 * \author Thomas M. Evans
 * \date   Tue Jul 18 16:07:35 2000
 * \brief  AMR_Layout class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_AMR_Layout_hh__
#define __mc_AMR_Layout_hh__

#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <iostream>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class AMR_Layout

 * \brief Provides cell connectivity information for 1-N meshes.
 
 * Class that defines the cell-face-cell connectivity for 1-to-N cell meshes.
 * 1-to-N describes meshes whose cells may have up to N cell neighbors across
 * a face.  AMR_Layout, like Layout, is independent of the geometry of the
 * mesh.

 * This class is designed for use in particular types of meshes; particularly
 * those that have 1-to-N cells across faces.  Therefore, we have not
 * designed this into a inheritance hierarchy with Layout.  Each mesh will
 * use the appropriate services that are required by a particular Layout.
 * The actual interface to the outside world is through the mesh type, not
 * the layout class.  The only classes that see the layout are the mesh and
 * its builder; thus, there is no need to restrain the interface to a common
 * set of functions declared by a base class.

 * The ordering of cells across a face is row-wise starting from the lowest
 * coordinate: thus across the hi x face in XY, 1=lo y and 2=hi y.
 * Similarly, across the hi y face in XYZ 1=lo x lo z, 2=hi x lo z, 3=lo x hi
 * z, and 4=hi x hi z.  There is no formal requirement to follow this
 * methodology; however, for consistency one should attempt to use this
 * convention.

 * The AMR_Layout uses the convention of coarse faces.  For example, in a 2-D
 * RZ mesh each cell has 4 coarse faces.  However, if one of the sides of a
 * cell borders a refined cell then there are 2 fine faces (2 cell neighbors)
 * on that coarse face.  This presents a "hierarchical" notion of a refined
 * mesh.

 * As in Layout, the data is built into AMR_Layout through use of the
 * overloaded operator().

 * \sa Layout class for information on 1-to-1 layouts and mesh layouts in
 * general.

 */
// revision history:
// -----------------
// 0) original
// 1) 19-APR-01 : added packing capability
// 
//===========================================================================//

class AMR_Layout 
{
  public:
    // Forward declaration of pack class.
    struct Pack;

  public:
    // Important typedefs.
    typedef std::vector<int>               sf_int;
    typedef std::vector<sf_int>            vf_int;
    typedef std::vector<vf_int>            tf_int;
    typedef rtt_dsxx::SP<AMR_Layout>       SP_Layout;
    typedef rtt_dsxx::SP<AMR_Layout::Pack> SP_Pack;

  private:
    //! Cell-face data, indexing is 0:N-1.
    tf_int cell_face;

  public:
    // Default constructor.
    explicit AMR_Layout(int = 0);

    // Constructor given the dimension.
    AMR_Layout(int, int);

    // >>> ACCESSORS

    //! Get the number of cells in the layout.
    int num_cells() const { return cell_face.size(); }

    // Get the number of coarse faces on a cell.
    inline int num_faces(int) const;

    // Number of cells across a coarse face.
    inline int num_cells_across(int, int) const;

    // >>> MODIFIERS

    //! Set size of layout to number of cells.
    void set_size(int num_cells) { cell_face.resize(num_cells); }

    // Set the number of coarse faces for a given cell.
    inline void set_size(int, int);

    // Set the number of fine faces across a coarse face in a given cell.
    inline void set_size(int, int, int);

    // >>> OVERLOADED OPERATORS
    
    // Overloaded subscripting operators for assignment and retrieval.
    inline int    operator()(int, int, int) const;
    inline int&   operator()(int, int, int);
    inline sf_int operator()(int, int) const;

    // Pack function.
    SP_Pack pack(const sf_int &) const;

    // Overloaded equality operators.
    bool operator==(const AMR_Layout &) const;
    bool operator!=(const AMR_Layout &rhs) const { return !(*this == rhs); }

    // Diagnostic funcitions.
    void print(std::ostream &, int) const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// Overload operator for stream output.

std::ostream& operator<<(std::ostream &, const AMR_Layout &);

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR AMR_LAYOUT
//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of coarse faces on a cell.
 */
int AMR_Layout::num_faces(int cell) const
{
    Require (cell > 0 && cell <= cell_face.size());
    return cell_face[cell-1].size();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of cells across a coarse face.
 */
int AMR_Layout::num_cells_across(int cell, int coarse_face) const
{
    Require (cell > 0 && cell <= cell_face.size());
    Require (coarse_face > 0 && coarse_face <= cell_face[cell-1].size());
    return cell_face[cell-1][coarse_face-1].size();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the number of coarse faces for a given cell.

 * By default this function gives 1 cell across each coarse face.

 * \param cell cell index
 * \param num_coarse_faces number of coarse faces for cell
 */
void AMR_Layout::set_size(int cell, int num_coarse_faces)
{
    cell_face[cell-1].resize(num_coarse_faces, sf_int(1));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the number of cells across a coarse face.
 
 * This allows the user to set the number of cells (fine faces) across a
 coarse face.

 * \param cell cell index
 * \param face coarse face index
 * \param num_cells_across number of cells across coarse face

 */
void AMR_Layout::set_size(int cell, int face, int num_cells_across)
{
    cell_face[cell-1][face-1].resize(num_cells_across);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access the cell across a coarse/fine face combination.
 * \param cell cell index
 * \param coarse_face coarse face index
 * \param fine_face fine face index
 * \return cell index on other side of cell, coarse_face, fine_face triplet
 */
int AMR_Layout::operator()(int cell, int coarse_face, int fine_face) const
{
    Require (cell > 0 && cell <= cell_face.size());
    Require (coarse_face > 0 && coarse_face <= cell_face[cell-1].size());
    Require (fine_face > 0 &&
	     fine_face <= cell_face[cell-1][coarse_face-1].size());
    return cell_face[cell-1][coarse_face-1][fine_face-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mutable access for the cell across a coarse/fine face combination.

 * This overloaded operator allows the user to change the value of the cell
 * across a face.  This operator should be used to build the AMR_Layout.

 * \param cell cell index
 * \param coarse_face coarse face index
 * \param fine_face fine face index

 * \return mutable cell index on other side of cell, coarse_face, fine_face
 * triplet

 */
int& AMR_Layout::operator()(int cell, int coarse_face, int fine_face)
{
    Require (cell > 0 && cell <= cell_face.size());
    Require (coarse_face > 0 && coarse_face <= cell_face[cell-1].size());
    Require (fine_face > 0 &&
	     fine_face <= cell_face[cell-1][coarse_face-1].size());
    return cell_face[cell-1][coarse_face-1][fine_face-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a vector of cells across a coarse face.

 * This is a non-mutable look-up option for determining the cells across a
 * coarse face.

 * \param cell cell index
 * \param coarse_face coarse face index
 * \return vector of cells across the coarse face

 */
AMR_Layout::sf_int AMR_Layout::operator()(int cell, int coarse_face) const
{
    Require (cell > 0 && cell <= cell_face.size());
    Require (coarse_face > 0 && coarse_face <= cell_face[cell-1].size());
    return cell_face[cell-1][coarse_face-1];
}

//===========================================================================//
/*!  
 * \struct Layout::Pack 
 
 * \brief Nested class for packing layouts into raw data for writing
 * or communication.
 
 */
//===========================================================================//

struct AMR_Layout::Pack
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

#endif                          // __mc_AMR_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of mc/AMR_Layout.hh
//---------------------------------------------------------------------------//
