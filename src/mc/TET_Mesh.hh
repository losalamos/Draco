//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Mesh.hh
 * \author H. Grady Hughes
 * \date   Tue Jan 18 10:15:43 MST 2000
 * \brief  Tetrahedral mesh class header file.
 */
//---------------------------------------------------------------------------//

#ifndef __mc_TET_Mesh_hh__
#define __mc_TET_Mesh_hh__

#include "ThreeVector.hh"
#include "Coord_sys.hh"
#include "Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace rtt_mc
{

// stl namespaces
using std::fill;
using std::min_element;
using std::ostream;
using std::pow;
using std::endl;
using std::string;

// draco namespaces
using rtt_rng::Sprng;
using dsxx::SP;

//___________________________________________________________________________//
/*!
 * \class TET_Mesh
 *
 * \brief General singly-connected tetrahedral mesh.
 *
 * Among private data, "internal list" and "internal numbers" refer to C-style
 * numbering from zero. External lists and numbers are FORTRAN-style, starting
 * from one.  In general, the arguments and return values of public functions
 * are external numbers.  Private data and implementation-detail calculations
 * involve internal numbers.
 *
 * A tetrahedral mesh is an inherently 3-dimensional Cartesian object.
 * The TET_Mesh class relies explicitly on this specialization at various
 * points in the code.
 */
// revision history:
// -----------------
//  0)   Original : Committed 2000_01_27.
//  1) 2000_01_28 : Added public functions to TET_Mesh:
//                          bool   in_cell(const SF_DOUBLE &, int);
//                          int    get_cell(const SF_DOUBLE &);
//                          double get_min_db(const SF_DOUBLE &, int);
//
//___________________________________________________________________________//

class TET_Mesh
{
 private:

    //! Private copy constructor: can't copy or assign a mesh.
    TET_Mesh(const TET_Mesh &);

    //! Private assignment operator: can't copy or assign a mesh.
    TET_Mesh& operator=(const TET_Mesh &);

    //! Typedef for scalar field of integers.
    typedef std::vector<int> SF_INT;

    //! Typedef for scalar field of doubles.
    typedef std::vector<double> SF_DOUBLE;

    //! Typedef for scalar field of ThreeVectors.
    typedef std::vector<ThreeVector> SF_THREEVECTOR;

    //! Typedef for vector field of integers.
    typedef std::vector< std::vector<int> > VF_INT;

    //! Typedef for vector field of doubles.
    typedef std::vector< std::vector<double> > VF_DOUBLE;

    //! The TET_Mesh is inherently 3-dimensional.
    static const int THREE = 3;

    //! A TET_Mesh cell always has 4 faces and 4 vertices.
    static const int FOUR = 4;

    //___________________________________________//
    // Beginning of private data of class TET_Mesh.

    //! Base class reference to a derived coordinate system class.
    SP<Coord_sys> coord;

    // Layout of mesh.
    Layout layout;

    /*!
     * vertex_vector[vertex#] == ThreeVector for the vertex labeled vertex#.
     *
     * vertex# == (0, 1, 2, 3, 4, ...) for an internal list of all vertices.
     */
    SF_THREEVECTOR vertex_vector;

    /*!
     * cells_vertices[cell#][cell_vertex#] == internal numbers of the four
     * vertices belonging to the given cell, labeled cell#.
     *
     * cell# == (0, 1, 2, 3, 4, ...) for an internal list of all cells.
     *
     * cell_vertex# == (0, 1, 2, 3) for each tetrahedral cell.
     *
     * When faces are identified, it is understood that face #0 is opposite
     * vertex #0, and is bounded by vertices #1, #2, #3, with corresponding
     * definitions for faces #1, #2, and #3.
     *
     * When cells_vertices[][] is constructed, it will be arranged that the
     * vertices labeled (0, 1, 2, 3) are in the proper order for predictable
     * inward- and outward-normal directions.
     */
    VF_INT cells_vertices;

    //! Flag to indicate whether this is a submesh.
    bool submesh;

    //_____________________________________//
    // End of private data of class TET_Mesh.

    void Valid(int cell) const
    { Require ( cell >= 1 && cell <= cells_vertices.size() ); }

    void Valid(int cell, int face) const
    {
        Require ( cell >= 1 && cell <= cells_vertices.size() );
        Require ( face >= 1 && face <= FOUR );
    }

    const ThreeVector get_outward_cross(int, int) const;

    const ThreeVector get_inward_cross(int, int) const;

 public:

    //! TET_Mesh constructor.
    TET_Mesh(SP<Coord_sys>,Layout &,SF_THREEVECTOR &,VF_INT &,bool = false);

    //! Forward declaration of cell-centered scalar fields.
    template<class T> class CCSF;

    //! Forward declaration of cell-centered vector fields.
    template<class T> class CCVF;

    //! Return the number of cells in the mesh.
    int num_cells() const { return layout.num_cells(); }

    // Determine whether a given position is inside a given cell.
    bool in_cell(const SF_DOUBLE &, int) const;

    //___________________________________________________//
    // Services required by all mesh types used in JAYENNE.

    int next_cell(int cell, int face) const
    { Valid(cell, face); return layout(cell, face); }

    double get_db(const SF_DOUBLE &, const SF_DOUBLE &, int, int &) const;

    const SF_DOUBLE get_normal(int, int) const;

    const SF_DOUBLE get_normal_in(int, int) const;

    double volume(int) const;

    double face_area(int, int) const;

    const SF_INT get_surcells(string) const;

    int get_bndface(string, int) const;

    const VF_DOUBLE get_vertices(int) const;

    const VF_DOUBLE get_vertices(int, int) const;

    const SF_DOUBLE sample_pos(int, Sprng &) const;

    const SF_DOUBLE sample_pos(int, Sprng &, SF_DOUBLE,
                     double) const;

    const SF_DOUBLE sample_pos_on_face(int, int, Sprng &) const;

    //! Return a vector<int> list of a cell's neighbors.
    const SF_INT get_neighbors(int cell) const
    {
        Valid(cell);
        SF_INT neighbors(layout.num_faces(cell));
        for (int face = 1; face <= neighbors.size(); face++)
        neighbors[face-1] = layout(cell, face);
        return neighbors;
    }

    bool full_Mesh() const { return !submesh; }

    // Overloaded comparison operators.
    bool operator==(const TET_Mesh &) const;
    bool operator!=(const TET_Mesh &rhs) const { return !(*this == rhs); }

    //_____________________________________//
    // End of required JAYENNE mesh services.

    // Find the cell comtaining a given position.
    int get_cell(const SF_DOUBLE &) const;

    // Get minimum distance to cell boundary, regardless of direction.
    double get_min_db(const SF_DOUBLE &, int) const;

    // References to imbedded objects and data required for Parallel_Building.
    // More may be added later.
    const Layout&    get_Layout() const  { return layout; }
    const Coord_sys& get_Coord() const   { return *coord; }
    SP<Coord_sys>    get_SPCoord() const { return coord; }

};  // end class TET_Mesh

//___________________________________________________________________________//
/*!
 * \class TET_Mesh::CCSF
 *
 * \brief Cell-centered scalar fields: one object per cell.
 *
 * This class does not require copy constructors or assignment operators
 * as the SP<> and vector<> classes can do assignment.
 *
 * The cell centered fields are used only outside of the MESH-type class.
 *
 * Empty fields cannot be built (the mesh must have at least one cell).
 */
template<class T>
class TET_Mesh::CCSF
{
 private:

    // Typedefs for notational convenience.
//  typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
//  typedef typename std::vector<T>::pointer pointer;
//  typedef typename std::vector<T>::const_pointer const_pointer;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::size_type size_type;

    //! Smart pointer back to underlying TET_Mesh.
    SP<TET_Mesh> mesh;

    //! Data in scalar field  --  one object per cell.
    vector<T> data;

 public:

    // Inline explicit constructor.
    inline explicit CCSF(SP<TET_Mesh>);

    // Constructor for automatic initialization.
    inline CCSF(SP<TET_Mesh>, const vector<T> &);

    //! Return reference to mesh.
    const TET_Mesh& get_Mesh() const { return *mesh; }

    /*!
     * Subscripting operations.
     *
     * We use () to indicate external cell numbers (1, 2, 3, 4, 5, etc.),
     * and [] to indicate internal identifiers (0, 1, 2, 3, 4, etc.).
     */

    //! \brief Constant overloaded () operator.
    const_reference operator()(int cell) const { return data[cell-1]; }

    //! \brief Assignment overloaded () operator.
    reference operator()(int cell) { return data[cell-1]; }

    //! \brief Constant overloaded [] operator.
    const_reference operator[](int cell_) const { return data[cell_]; }

    //! \brief Assignment overloaded [] operator.
    reference operator[](int cell_) { return data[cell_]; }

    // STL style functions.

    //! Iterators from first cell to one-past-last cell.
    iterator begin() { return data.begin(); }

    const_iterator begin() const { return data.begin(); }

    iterator end() { return data.end(); }

    const_iterator end() const { return data.end(); }

    //! \brief For scalar fields, return the number of cells.
    size_type size() const { return data.size(); }

    bool empty() const { return data.empty(); }

};  // end class TET_Mesh::CCSF

//---------------------------------------------------------------------------//
//  TET_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//

//___________________________________________________________________________//
/*!
 * \brief       CCSF explicit constructor.
 * \param mesh_ Smart pointer to the underlying mesh for this scalar field.
 */
template<class T>
inline TET_Mesh::CCSF<T>::CCSF(SP<TET_Mesh> mesh_)
    : mesh(mesh_), data(mesh->num_cells())
{
    Require (mesh);
    Ensure (!empty());
}

//___________________________________________________________________________//
/*!
 * \brief       CCSF constructor copying an existing scalar field.
 * \param mesh_ Smart pointer to the underlying mesh for the new scalar field.
 * \param array The existing scalar field (presumably) on the same mesh.
 *
 * Constructor for automatic initialization.
 */
template<class T>
inline TET_Mesh::CCSF<T>::CCSF(SP<TET_Mesh> mesh_, const vector<T> &array)
    : mesh(mesh_), data(array)
{
    Require (mesh);
    Ensure  (data.size() == mesh->num_cells());
    Ensure  (!empty());
}

//---------------------------------------------------------------------------//
//  end of TET_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//

//___________________________________________________________________________//
/*!
 * \class TET_Mesh::CCVF
 *
 * \brief Cell-centered vector fields.
 *
 * This class does not require copy constructors or assignment operators
 * as the SP<> and vector<> classes can do assignment.
 *
 * The cell centered fields are used only outside of the MESH-type class.
 *
 * Empty fields cannot be built (the mesh must have at least one cell).
 */
template<class T>
class TET_Mesh::CCVF
{
 private:

    // Typedefs for notational convenience.
//  typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
//  typedef typename std::vector<T>::pointer pointer;
//  typedef typename std::vector<T>::const_pointer const_pointer;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::size_type size_type;

    //! Smart pointer back to underlying TET_Mesh.
    SP<TET_Mesh> mesh;

    //! Data in vector field  --  one object per (dimension, cell) pairing.
    vector< vector<T> > data;

 public:

    // Inline explicit constructor.
    inline explicit CCVF(SP<TET_Mesh>);

    // Constructor for automatic initialization.
    inline CCVF(SP<TET_Mesh>, const vector<vector<T> > &);

    //! Return reference to mesh.
    const TET_Mesh& get_Mesh() const { return *mesh; }

    /*!
     * Subscripting operations.
     *
     * We use () to indicate external cell numbers (1, 2, 3, 4, 5, etc.),
     * and external dimensions (1, 2, 3 for X, Y, Z).
     */

    //! \brief Constant overloaded () operator.
    const_reference operator()(int dim, int cell) const
    { return data[dim-1][cell-1]; }

    //! \brief Assignment overloaded () operator.
    reference operator()(int dim, int cell)
    { return data[dim-1][cell-1]; }

    // Get a vector associated with a cell.
    inline vector<T> operator()(int) const;

    // STL style functions.

    //! Iterators without argument refer to dimensions (X,Y,Z) - NOT cells.
    const_iterator begin() const { return data.begin(); }

    iterator begin() { return data.begin(); }

    const_iterator end() const { return data.end(); }

    iterator end() { return data.end(); }

    //! Iterators with arguments allow iteration over cells.
    const_iterator begin(int i) const
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].begin(); }

    iterator begin(int i)
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].begin(); }

    const_iterator end(int i) const
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].end(); }

    iterator end(int i)
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].end(); }

    //! Without argument, size() should always return 3.
    size_type size() const { return data.size(); }

    //! size(i) returns the number of entries (cells) in the i-direction.
    size_type size(int i) const
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].size(); }

    //! Without argument, empty() checks the entire field.
    bool empty() const { return data.empty(); }

    //! empty(i) checks absence of entries in the i-dimension of the field.
    bool empty(int i) const
    { Require ( i >= 1 && i <= data.size() ); return data[i-1].empty(); }

};  // end class TET_Mesh::CCVF

//---------------------------------------------------------------------------//
//  TET_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//

//___________________________________________________________________________//
/*!
 * \brief       CCVF explicit constructor.
 * \param mesh_ Smart pointer to the underlying mesh for this vector field.
 */
template<class T>
inline TET_Mesh::CCVF<T>::CCVF(SP<TET_Mesh> mesh_)
    : mesh(mesh_), data(mesh->get_Coord().get_dim())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i < mesh->get_Coord().get_dim(); i++)
        data[i].resize(mesh->num_cells());
}

//___________________________________________________________________________//
/*!
 * \brief       CCVF constructor copying an existing vector field.
 * \param mesh_ Smart pointer to the underlying mesh for the new vector field.
 * \param array The existing vector field (presumably) on the same mesh.
 *
 * Constructor for automatic initialization.
 */
template<class T>
inline TET_Mesh::CCVF<T>::CCVF(SP<TET_Mesh> mesh_,
                  const vector<vector<T> > &array)
    : mesh(mesh_), data(array)
{
    Ensure (data.size() == mesh->get_Coord().get_dim());

    for (int dim = 0; dim < mesh->get_Coord().get_dim(); dim++)
        Ensure (data[dim].size() == mesh->num_cells());
}

//___________________________________________________________________________//
/*!
 * \brief      Overload () operator to get a vector associated with a cell.
 * \param cell External number of the given cell.
 * \return     The vector of values represented by the first index of the CCVF.
 *
 * For TET_Mesh, the returned vector is 3-dimensional.
 */
template<class T>
inline vector<T> TET_Mesh::CCVF<T>::operator()(int cell) const
{
    // declare return vector
    vector<T> x;

    // loop through dimensions and make return vector for this cell
    for (int i = 0; i < data.size(); i++)
        x.push_back(data[i][cell-1]);

    Ensure (x.size() == data.size());
    return x;
}

//---------------------------------------------------------------------------//
//  end of TET_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//

}   // end namespace rtt_mc

#endif  // __mc_TET_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of mc/TET_Mesh.hh
//---------------------------------------------------------------------------//
