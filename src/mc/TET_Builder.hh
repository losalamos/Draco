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

namespace rtt_mc
{

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
//  0)   Original : Committed 2000-02-08
//  1) 2000-02-12 : Completed namespace issues, bug fixes, and elimination
//                  of using declarations.
//  2) 2000-04-10 : Rewritten to be consistent with new meshReader classes.
//
//___________________________________________________________________________//

class TET_Builder
{
 private:

    //! Typedef for scalar field of integers.
    typedef std::vector<int> SF_INT;

    //! Typedef for scalar field of ThreeVectors.
    typedef std::vector<ThreeVector> SF_THREEVECTOR;

    //! Typedef for vector field of integers.
    typedef std::vector< std::vector<int> > VF_INT;

    //! Typedef for vector field of doubles.
    typedef std::vector< std::vector<double> > VF_DOUBLE;

    //! The TET_Mesh is inherently 3-dimensional and its faces have 3 vertices.
    static const int THREE = 3;

    //! A TET_Mesh cell always has 4 faces and 4 vertices.
    static const int FOUR = 4;

    //______________________________________________//
    // Beginning of private data of class TET_Builder.
    // Vector indices are all internal numbers.

    //! Coordinate system identifying string.
    std::string coord_system;

    //! Collection of vertices: nodes_coords[v][dim]==dim-coordinate of node v.
    VF_DOUBLE nodes_coords;

    /*!
     * Parent nodes for each vertex: parent[v] == external number of parent
     * node for the internally-numbered node v.
     * These will be used only to establish connectivity (layout). They need
     * not be converted to internal numbers, as long as care is taken with
     * their consistent use.
     */
    SF_INT parent;

    /*!
     * Nodes bounding each cell: cells_vertices[c][0..3] == internal # of node.
     * These can be provided in any order.  The TET_Builder constructor will
     * adjust them to provide the right-handed sense needed by TET_Mesh.  Also
     * the interface (e.g. file reader) is expected to provide internal numbers
     * of the bounding nodes, so no adjustment is needed on construction.
     */
    VF_INT cells_vertices;

    //! Flag to indicate whether this is a submesh.
    bool submesh;

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
    nodes_coords   = interface->get_node_coords();

    // TET_Mesh objects are 3-D.
    for (int node_ = 0 ; node_ < nodes_coords.size() ; node_++)
        Check (nodes_coords[node_].size() == THREE);

    cells_vertices = interface->get_element_nodes();

    VF_INT::iterator first_cell = cells_vertices.begin();
    for ( ; first_cell != cells_vertices.end() ; first_cell++)
        if ((*first_cell).size() != THREE)
            break;

    Require (first_cell != cells_vertices.begin());
    Require (first_cell != cells_vertices.end());

    cells_vertices.erase(cells_vertices.begin(),first_cell);

    // TET_Mesh objects have four vertices.
    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
        Check (cells_vertices[cell_].size() == FOUR);

    parent.resize(nodes_coords.size());
    for (int node_ = 0 ; node_ < nodes_coords.size() ; node_++)
        parent[node_] = node_ + 1;

    // To be properly constructed, cells_vertices[][] must contain the
    // internal numbers of the four bounding vertices of the cells.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
        for (int ver_ = 0 ; ver_ < FOUR ; ver_++)
        {
            Ensure ( cells_vertices[cell_][ver_] >= 0 );
            Ensure ( cells_vertices[cell_][ver_] < nodes_coords.size() );
        }

    // To be properly constructed, cells_vertices[][] must be in a predefined
    // order, so that later sense-of-surface operations will be predictable.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
    {
        int v0 = cells_vertices[cell_][0];
        int v1 = cells_vertices[cell_][1];
        int v2 = cells_vertices[cell_][2];
        int v3 = cells_vertices[cell_][3];

        ThreeVector R0(nodes_coords[v0][0],nodes_coords[v0][1],
                       nodes_coords[v0][2]);
        ThreeVector R1(nodes_coords[v1][0],nodes_coords[v1][1],
                       nodes_coords[v1][2]);
        ThreeVector R2(nodes_coords[v2][0],nodes_coords[v2][1],
                       nodes_coords[v2][2]);
        ThreeVector R3(nodes_coords[v3][0],nodes_coords[v3][1],
                       nodes_coords[v3][2]);

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
