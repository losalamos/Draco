//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Mesh.hh
 * \author H. Grady Hughes
 * \date   Tue Jan 18 10:15:43 MST 2000
 * \brief  Tetrahedral mesh class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_TET_Mesh_hh__
#define __mc_TET_Mesh_hh__

#include "ThreeVector.hh"
#include "XYZCoord_sys.hh"
#include "Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "viz/Ensight_Translator.hh"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <set>
#include <map>

namespace rtt_mc
{

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
//  0)   Original: Committed 2000-01-27.
//  1) 2000-01-28: Added public functions to TET_Mesh:
//                         bool   in_cell(const sf_double &, int);
//                         int    get_cell(const sf_double &);
//                         double get_min_db(const sf_double &, int);
//  2) 2000-01-30: Added subscripting, iterators, size, and empty functions
//                 to class TET_Mesh::CCVF.
//  3) 2000-02-12: Added TET_Mesh member functions get_cell_types() and
//                 get_point_coord() for Ensight support.  Also completed all
//                 namespace issues and elimination of using declarations.
//  4) 2000-04-25: Renamed in_cell() to in_open_cell() and added the related
//                 function in_closed_cell().
//  5) 2000-04-26: Added private data vf_int sides_vertices, and dealt with
//                 sides_vertices in the public constructor.
//  6) 2000-05-03: TET_Builder, TET_Mesh, and their test files now use the
//                 get_node_coord_units(), get_node_sets(), get_element_sets(),
//                 and get_title() services of the Mesh_Reader base class.
//                 At the top level (TET_Mesh), the get_element_sets() services
//                 will later be replaced by side- and cell-specific data
//                 structures.
//  7) 2000-06-08: Information from the interface service get_element_sets()
//                 is now converted to two separate maps, side_sets and
//                 cell_sets, and used to initialize data members of the
//                 TET_Mesh class.  The TET_Mesh class no longer has knowledge
//                 of element_sets.  New diagnostic functions print_node_sets,
//                 print_side_sets, and print_cell_sets are added to TET_Mesh.
//  8) 2000-06-19: Corrected a major bug in get_db() and get_min_db().
//  9) 2000-07-23: Initial implementation of sample_pos() for linearly
//                 interpolated temperature**4 given the temperatures**4 on the
//                 vertices of the cells.  Its interface will probably need to
//                 be changed; and it is less efficient than it might be,
//                 because of recalculation of certain quantities that we will
//                 probably decide to store.
// 10) 2000-11-21: Implemented print_vertex_vector() member function, which
//                 was the last part of the overall print_mesh() function.
// 11) 2000-11-30: Added get_cell_pair(), the third graphic service for the
//                 TET_Mesh class, to support Ensight_Translator.
//
//___________________________________________________________________________//

class TET_Mesh
{
 public:

    //! Forward declaration of pack class.
    struct Pack;

    //! Forward declaration of cell-centered scalar fields.
    template<class T> class CCSF;

    //! Forward declaration of cell-centered vector fields.
    template<class T> class CCVF;

    // Public-interface typedefs.
    typedef rtt_dsxx::SP<TET_Mesh>               SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>              SP_Coord_sys;
    typedef rtt_dsxx::SP<TET_Mesh::Pack>         SP_Pack;
    typedef std::string                          std_string;

    // Public-interface typedefs for fields of standard types.
    typedef std::vector<int>                     sf_int;
    typedef std::vector< std::vector<int> >      vf_int;
    typedef std::vector<double>                  sf_double;
    typedef std::vector< std::vector<double> >   vf_double;
    typedef std::vector<std_string>              sf_string;

    // Typedefs to cell-centered fields.
    typedef CCSF<double>                         CCSF_double;
    typedef CCSF<int>                            CCSF_int;
    typedef CCVF<double>                         CCVF_double;
    typedef CCVF<int>                            CCVF_int;
    typedef CCSF<std_string>                     CCSF_string;

 private:

    //! Private copy constructor: can't copy or assign a mesh.
    TET_Mesh(const TET_Mesh &);

    //! Private assignment operator: can't copy or assign a mesh.
    TET_Mesh& operator=(const TET_Mesh &);

    //! Typedef for scalar field of ThreeVectors.
    typedef std::vector<ThreeVector> sf_ThreeVector;

    //! Typedef for a standard set of integers.
    typedef std::set<int> SetInt;

    //! Typedef for a map linking strings to sets of integers.
    typedef std::map< std::string, SetInt > MAP_String_SetInt;

    //! The TET_Mesh is inherently 3-dimensional and its faces have 3 vertices.
    static const int THREE = 3;

    //! A TET_Mesh cell always has 4 faces and 4 vertices.
    static const int FOUR = 4;

    //___________________________________________//
    // Beginning of private data of class TET_Mesh.

    //! Mesh title.
    std::string title;

    //! Base class reference to a derived coordinate system class.
    SP_Coord_sys coord;

    // Layout of mesh.
    Layout layout;

    /*!
     * vertex_vector[vertex#] == ThreeVector for the vertex labeled vertex#.
     *
     * vertex# == (0, 1, 2, 3, 4, ...) for an internal list of all vertices.
     */
    sf_ThreeVector vertex_vector;

    //! Coordinate system units (e.g. "cm").
    std::string node_coord_units;

    //! Associate sets of nodes with characteristics identified by strings.
    MAP_String_SetInt node_sets;

    //! Associate sets of sides with characteristics identified by strings.
    MAP_String_SetInt side_sets;

    //! Associate sets of cells with characteristics identified by strings.
    MAP_String_SetInt cell_sets;

    /*!
     * sides_vertices[side#][side_vertex#] == internal numbers of the three
     * vertices belonging to the given side, labeled side#.
     *
     * side# == (0, 1, 2, 3, 4, ...) for an internal list of all sides.
     *
     * side_vertex# == (0, 1, 2) for each triangular side.
     *
     * No assumptions are currently made about the relation of any given side
     * to any particular cell face, nor about the ordering of the vertices of
     * a side.
     */
    vf_int sides_vertices;

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
     * Before cells_vertices[][] is constructed, it will have been arranged
     * that the vertices labeled (0, 1, 2, 3) are in the proper order for
     * predictable inward- and outward-normal directions.
     */
    vf_int cells_vertices;

    //! Flag to indicate whether this is a submesh.
    bool submesh;

    //_____________________________________//
    // End of private data of class TET_Mesh.

    //! \brief Make sure that an external cell number is valid.
    void Valid(int cell) const
    { Require ( cell >= 1 && cell <= cells_vertices.size() ); }

    //! \brief Make sure that external cell and face numbers are valid.
    void Valid(int cell, int face) const
    {
        Require ( cell >= 1 && cell <= cells_vertices.size() );
        Require ( face >= 1 && face <= FOUR );
    }

    const ThreeVector get_outward_cross(int, int) const;

    const ThreeVector get_inward_cross(int, int) const;

    const sf_double get_barycentric_coords(const sf_double &, int) const;

 public:

    //! TET_Mesh constructor.
    TET_Mesh(std::string &, SP_Coord_sys, Layout &,
        sf_ThreeVector &, std::string &, MAP_String_SetInt &,
        MAP_String_SetInt &, MAP_String_SetInt &, vf_int &, vf_int &,
        bool = false);

    //! \brief Return the number of cells in the mesh.
    int num_cells() const { return layout.num_cells(); }

    // Determine whether a given position is inside a given cell.
    bool in_open_cell(const sf_double &, int) const;

    // Similarly, inside a cell or on its boundary.
    bool in_closed_cell(const sf_double &, int) const;

    //____________________________________________________//
    // Services required by all mesh types used in JAYENNE.

    //! \brief Cell on the other side of "face" from "cell"
    int next_cell(int cell, int face) const
    { Valid(cell, face); return layout(cell, face); }

    double get_db(const sf_double &, const sf_double &, int, int &) const;

    const sf_double get_normal(int, int) const;

    const sf_double get_normal_in(int, int) const;

    double volume(int) const;

    double face_area(int, int) const;

    const sf_int get_surcells(std::string) const;

    int get_bndface(std::string, int) const;

    const vf_double get_vertices(int) const;

    const vf_double get_vertices(int, int) const;

    const sf_double sample_pos(int, rtt_rng::Sprng &) const;

    const sf_double sample_pos(int, rtt_rng::Sprng &, const sf_double &) const;

    const sf_double sample_pos_on_face(int, int, rtt_rng::Sprng &) const;

    //! \brief Return a vector<int> list of a cell's neighbors.
    const sf_int get_neighbors(int cell) const
    {
        Valid(cell);
        sf_int neighbors(layout.num_faces(cell));
        for (int face = 1; face <= neighbors.size(); face++)
        neighbors[face-1] = layout(cell, face);
        return neighbors;
    }

    bool full_Mesh() const { return !submesh; }

    // Overloaded comparison operators.
    bool operator==(const TET_Mesh &) const;
    bool operator!=(const TET_Mesh &rhs) const { return !(*this == rhs); }

    //______________________________________//
    // End of required JAYENNE mesh services.

    //_____________________________________//
    // Services required for graphics dumps.

    //! \brief Return the cell type for each cell in the mesh.
    sf_int get_cell_types() const
    {
        sf_int cell_type(layout.num_cells());
        std::fill(cell_type.begin(), cell_type.end(),
            rtt_viz::four_node_tetrahedron);

        return cell_type;
    }

    /*!
     * \brief Return vector coordinates for each vertex in the mesh.
     *
     * return_coords[vertex#][0,1,2] == (X,Y,Z) coordinate of vertex#.
     * vertex# == (0, 1, 2, 3, 4, ...) for an internal list of all vertices.
     */
    vf_double get_point_coord() const
    {
        vf_double return_coords(vertex_vector.size());
        for (int v_ = 0 ; v_ < vertex_vector.size() ; v_++)
        {
            return_coords[v_].push_back(vertex_vector[v_].get_x());
            return_coords[v_].push_back(vertex_vector[v_].get_y());
            return_coords[v_].push_back(vertex_vector[v_].get_z());
        }
        return return_coords;
    }

    /*!
     * \brief Return external vertex numbers for each vertex of each cell.
     *
     * cell_pair[cell#][cell_vertex#] == external number of the given
     * vertex belonging to the given cell.
     *
     * cell# == (0, 1, 2, 3, 4, ...) for an internal list of all cells.
     * cell_vertex# == (0, 1, 2, 3) for each tetrahedral cell.
     */
    vf_int get_cell_pair() const
    {
        vf_int cell_pair(cells_vertices);
        for (int i = 0; i < cell_pair.size(); ++i)
            for (int j = 0; j < cell_pair[i].size(); ++j)
                ++cell_pair[i][j];
        return cell_pair;
    }

    //______________________________//
    // End of graphics dump services.

    //____________________//
    // Diagnostic services.

    //! Print the mesh title.
    void print_title(std::ostream &) const;

    //! Print the layout.
    void print_layout(std::ostream &) const;

    //! Print the vertex_vector.
    void print_vertex_vector(std::ostream &) const;

    //! Print the node_coord_units.
    void print_node_coord_units(std::ostream &) const;

    //! Print the node_sets.
    void print_node_sets(std::ostream &) const;

    //! Print the side_sets.
    void print_side_sets(std::ostream &) const;

    //! Print the cell_sets.
    void print_cell_sets(std::ostream &) const;

    //! Print the sides_vertices.
    void print_sides_vertices(std::ostream &) const;

    //! Print the cells_vertices.
    void print_cells_vertices(std::ostream &) const;

    //! Print the submesh status.
    void print_submesh(std::ostream &) const;

    //! Print everything about the mesh.
    void print_mesh(std::ostream &) const;

    //___________________________//
    // End of diagnostic services.

    // Find the cell comtaining a given position.
    int get_cell(const sf_double &) const;

    // Get minimum distance to cell boundary, regardless of direction.
    double get_min_db(const sf_double &, int) const;

    // References to imbedded objects and data required for Parallel_Building.
    // More may be added later.
    const Layout&     get_Layout() const  { return layout; }
    const Coord_sys&  get_Coord() const   { return *coord; }
    SP_Coord_sys      get_SPCoord() const { return coord; }

    //! TET_Mesh member function to return a Pack object holding a TET_Mesh.
    SP_Pack pack() const;

};  // end class TET_Mesh

/*!
 * \struct TET_Mesh::Pack
 *
 * \brief Structure to hold a packed tetrahedral mesh.
 */
struct TET_Mesh::Pack
{
 private:
    // Disallow assignment.
    const Pack& operator=(const Pack &);

    //! TET_Mesh is inherently 3-dimensional; its faces have 3 vertices.
    static const int THREE = 3;

    //! A TET_Mesh cell always has 4 faces and 4 vertices.
    static const int FOUR = 4;

    //! Typedef for vector field of integers.
    typedef std::vector< std::vector<int> > vf_int;

    //! Typedef for a standard set of integers.
    typedef std::set<int> SetInt;

    //! Typedef for a map linking strings to sets of integers.
    typedef std::map< std::string, SetInt > MAP_String_SetInt;

    // Private data members.

    // Collected double data.
    int dsize;
    double *ddata;

    // Collected integer data.
    int isize;
    int *idata;

    // Collected character data.
    int csize;
    char *cdata;

 public:

    // Constructor.  Note that after construction, the Pack object
    // owns the pointed-to data, and assumes the responsibility of
    // deletion.
    Pack(int, double *, int, int *, int, char *);

    // Copy constructor.
    Pack(const Pack &);

    // Destructor.
    ~Pack();

    // Unpacker.  Recovers a TET_Mesh from the Pack object.
    TET_Mesh::SP_Mesh unpack() const;

    // Printer for diagnostics.
    void print_pack(std::ostream &) const;

};  // end struct TET_Mesh::Pack

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
    SP_Mesh mesh;

    //! Data in scalar field  --  one object per cell.
    std::vector<T> data;

 public:

    // Inline explicit constructor.
    inline explicit CCSF(SP_Mesh);

    // Constructor for automatic initialization.
    inline CCSF(SP_Mesh, const std::vector<T> &);

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
inline TET_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_)
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
inline TET_Mesh::CCSF<T>::CCSF(SP_Mesh mesh_,
    const std::vector<T> &array) : mesh(mesh_), data(array)
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
    typedef T& reference;
    typedef const T& const_reference;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::size_type size_type;

    //! Smart pointer back to underlying TET_Mesh.
    SP_Mesh mesh;

    //! Data in vector field  --  one object per (dimension, cell) pairing.
    std::vector< std::vector<T> > data;

 public:

    // Inline explicit constructor.
    inline explicit CCVF(SP_Mesh);

    // Constructor for automatic initialization.
    inline CCVF(SP_Mesh, const std::vector< std::vector<T> > &);

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
    inline std::vector<T> operator()(int) const;

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
inline TET_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_)
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
inline TET_Mesh::CCVF<T>::CCVF(SP_Mesh mesh_,
                  const std::vector< std::vector<T> > &array)
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
inline std::vector<T> TET_Mesh::CCVF<T>::operator()(int cell) const
{
    // declare return vector
    std::vector<T> x;

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
