//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Builder.hh
 * \author H. Grady Hughes
 * \date   Mon Jan 31 13:23:43 MST 2000
 * \brief  Tetrahedral mesh builder class header file.
 */
//---------------------------------------------------------------------------//

#ifndef __mc_TET_Builder_hh__
#define __mc_TET_Builder_hh__

#include "TET_Mesh.hh"
#include "XYZCoord_sys.hh"
#include "meshReaders/Element_Definition.hh"

namespace rtt_mc
{

using rtt_meshReaders::Element_Definition;

using std::cerr;  // FOR DEBUGGING.
using std::endl;  // FOR DEBUGGING.

//___________________________________________________________________________//
/*!
 * \class TET_Builder
 *
 * \brief Builds a TET_Mesh object.
 */
// revision history:
// -----------------
//  0)   Original: Committed 2000-02-08
//  1) 2000-02-12: Completed namespace issues, bug fixes, and elimination
//                 of using declarations.
//  2) 2000-04-10: Rewritten to be consistent with new meshReader classes.
//  3) 2000-04-26: Modified to construct VF_INT sides_vertices, new private
//                 TET_Mesh data, and to make use of get_element_nodes() and
//                 get_element_types() from the Mesh_Reader base class.
//  4) 2000-05-03: TET_Builder, TET_Mesh, and their test files now use the
//                 get_node_coord_units(), get_node_sets(), get_element_sets(),
//                 and get_title() services of the Mesh_Reader base class.
//                 At the top level (TET_Mesh), the get_..._sets() services
//                 will later be replaced by side- and cell-specific data
//                 structures.
//
//___________________________________________________________________________//

class TET_Builder
{
 private:

    //! Typedef for scalar field of integers.
    typedef std::vector<int> SF_INT;

    //! Typedef for scalar field of ThreeVectors.
    typedef std::vector<ThreeVector> SF_THREEVECTOR;

    //! Typedef for scalar field of Element_Types.
    typedef std::vector<Element_Definition::Element_Type> SF_TYPE;

    //! Typedef for vector field of integers.
    typedef std::vector< std::vector<int> > VF_INT;

    //! Typedef for vector field of doubles.
    typedef std::vector< std::vector<double> > VF_DOUBLE;

    //! Typedef for a standard set of integers.
    typedef std::set<int> SetInt;

    //! Typedef for a map linking strings to sets of integers.
    typedef std::map< std::string, SetInt > MAP_String_SetInt;

    //! The TET_Mesh is inherently 3-dimensional and its faces have 3 vertices.
    static const int THREE = 3;

    //! A TET_Mesh cell always has 4 faces and 4 vertices.
    static const int FOUR = 4;

    //______________________________________________//
    // Beginning of private data of class TET_Builder.
    // Vector indices are all internal numbers.

    //! Coordinate system identifying string.
    std::string coord_system;

    //! Flag to indicate whether this is a submesh.
    bool submesh;

    //! Collection of vertices: node_coords[v][dim]==dim-coordinate of node v.
    VF_DOUBLE node_coords;

    //! Coordinate system units (e.g. "cm").
    std::string node_coord_units;

    //! Associate sets of nodes with characteristics identified by strings.
    MAP_String_SetInt node_sets;

    //! Associate sets of elements with characteristics identified by strings.
    MAP_String_SetInt element_sets;

    //! Mesh title.
    std::string title;

    /*!
     * Parent nodes for each vertex: parent[v] == external number of parent
     * node for the internally-numbered node v.
     * These will be used only to establish connectivity (layout). They need
     * not be converted to internal numbers, as long as care is taken with
     * their consistent use.
     */
    SF_INT parent;

    /*!
     * Nodes bounding each side: sides_vertices[s][0..2] == internal # of node.
     * No assumptions are currently made about the relation of any given side
     * to any particular cell face, nor about the ordering of the vertices of
     * a side.
     * The interface (e.g. file reader) is expected to provide internal numbers
     * of the bounding nodes, so no adjustment is needed on construction.
     */
    VF_INT sides_vertices;

    /*!
     * Nodes bounding each cell: cells_vertices[c][0..3] == internal # of node.
     * These can be provided in any order.  The TET_Builder constructor will
     * adjust them to provide the right-handed sense needed by TET_Mesh.  Also
     * the interface (e.g. file reader) is expected to provide internal numbers
     * of the bounding nodes, so no adjustment is needed on construction.
     */
    VF_INT cells_vertices;

    //________________________________________//
    // End of private data of class TET_Builder.

 public:

    //! Explicit constructor, templated on Interface Type.
    template<class IT>
    explicit TET_Builder(rtt_dsxx::SP<IT>);

    //! Build a TET_Mesh, using TET_Builder's private data.
    rtt_dsxx::SP<TET_Mesh> build_Mesh();

};  // end class TET_Builder

//___________________________________________________________________________//
/*!
 * \brief           Non-inline, templated TET_Builder constructor.
 * \param interface Smart pointer to an instance of an Interface Type.
 *
 * This version does not require that the interface class return information
 * on parent nodes, the designator of the coordinate system, or the submesh
 * flag.  The first two of these items probably should use services of the
 * interface class.  However, for the time being, coord_system and submesh
 * are initialized, and parent is set to the obvious default.
 */
template<class IT>
TET_Builder::TET_Builder(rtt_dsxx::SP<IT> interface)
    : coord_system("xyz"), submesh(false)
{
    Require (interface);

    // Get data arrays from the interface for TET_Mesh objects.
    node_coords   = interface->get_node_coords();

    node_coord_units = interface->get_node_coord_units();

    VF_INT element_nodes = interface->get_element_nodes();

    SF_TYPE element_types = interface->get_element_types();

    node_sets = interface->get_node_sets();
// Begin: For initial debugging.
    for (MAP_String_SetInt::iterator flag = node_sets.begin() ; 
        flag != node_sets.end() ; flag++)
        {
            cerr << (*flag).first << endl;
            int nnode = (*flag).second.size();
            cerr << "  size = " << nnode << endl;
            if (nnode > 0)
                for (SetInt::iterator i = (*flag).second.begin() ;
                                       i != (*flag).second.end(); i++)
                    cerr << *i << endl;
        }
// End: For initial debugging.

    element_sets = interface->get_element_sets();
// Begin: For initial debugging.
    for (MAP_String_SetInt::iterator flag = element_sets.begin() ; 
        flag != element_sets.end() ; flag++)
        {
            cerr << (*flag).first << endl;
            int nelem = (*flag).second.size();
            cerr << "  size = " << nelem << endl;
            if (nelem > 0)
                for (SetInt::iterator i = (*flag).second.begin() ;
                                       i != (*flag).second.end(); i++)
                    cerr << *i << endl;
        }
// End: For initial debugging.

    title = interface->get_title();

    // TET_Mesh objects are 3-D.
    for (int node_ = 0 ; node_ < node_coords.size() ; node_++)
        Check (node_coords[node_].size() == THREE);

    Check (element_types.size() == element_nodes.size());

    parent.resize(node_coords.size());
    for (int node_ = 0 ; node_ < node_coords.size() ; node_++)
        parent[node_] = node_ + 1;

    int num_sides = 0;
    int num_cells = 0;

    for (int elem = 0 ; elem < element_types.size() ; elem++)
        if (element_types[elem] == Element_Definition::TRI_3)
            num_sides++;
        else if (element_types[elem] == Element_Definition::TETRA_4)
            num_cells++;
        else
            Check (element_types[elem] == Element_Definition::NODE ||
                   element_types[elem] == Element_Definition::BAR_2);

    sides_vertices.resize(num_sides);
    cells_vertices.resize(num_cells);

    for (int elem = 0, cell_ = 0, side_ = 0 ;
             elem < element_types.size() ; elem++)
        if (element_types[elem] == Element_Definition::TRI_3)
            {
                sides_vertices[side_] = element_nodes[elem];
                side_++;
            }
        else if (element_types[elem] == Element_Definition::TETRA_4)
            {
                cells_vertices[cell_] = element_nodes[elem];
                cell_++;
            }

    // To be properly constructed, sides_vertices[][] must contain the
    // internal numbers of the three bounding vertices of the sides.

    for (int side_ = 0 ; side_ < sides_vertices.size() ; side_++)
        for (int ver_ = 0 ; ver_ < THREE ; ver_++)
        {
            Ensure ( sides_vertices[side_][ver_] >= 0 );
            Ensure ( sides_vertices[side_][ver_] < node_coords.size() );
        }

    // To be properly constructed, cells_vertices[][] must contain the
    // internal numbers of the four bounding vertices of the cells.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
        for (int ver_ = 0 ; ver_ < FOUR ; ver_++)
        {
            Ensure ( cells_vertices[cell_][ver_] >= 0 );
            Ensure ( cells_vertices[cell_][ver_] < node_coords.size() );
        }

    // To be properly constructed, cells_vertices[][] must be in a predefined
    // order, so that later sense-of-surface operations will be predictable.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
    {
        int v0 = cells_vertices[cell_][0];
        int v1 = cells_vertices[cell_][1];
        int v2 = cells_vertices[cell_][2];
        int v3 = cells_vertices[cell_][3];

        ThreeVector R0(node_coords[v0][0],node_coords[v0][1],
                       node_coords[v0][2]);
        ThreeVector R1(node_coords[v1][0],node_coords[v1][1],
                       node_coords[v1][2]);
        ThreeVector R2(node_coords[v2][0],node_coords[v2][1],
                       node_coords[v2][2]);
        ThreeVector R3(node_coords[v3][0],node_coords[v3][1],
                       node_coords[v3][2]);

        double disc = (R1 - R0).dot((R2 - R1).cross(R3 - R2));

        Ensure (disc != 0.0);

        if (disc < 0.0)
        {
            int temp = cells_vertices[cell_][2];
            cells_vertices[cell_][2] = cells_vertices[cell_][3];
            cells_vertices[cell_][3] = temp;
        }
    }

}   // end TET_Builder::TET_Builder(rtt_dsxx::SP<IT>)

}   // end namespace rtt_mc

#endif  // __mc_TET_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/TET_Builder.hh
//---------------------------------------------------------------------------//
