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

// using std::cerr;  // FOR DEBUGGING.
// using std::endl;  // FOR DEBUGGING.

// draco namespaces
using rtt_rng::Sprng;
using rtt_dsxx::SP;

//___________________________________________________________________________//
/*!
 * \class TET_Builder
 *
 * \brief Builds a TET_Mesh object.
 */
// revision history:
// -----------------
//  0)   Original : Committed 
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

    //! Collection of vertices: nodes[dim][v] == dim-coordinate of node v.
    VF_DOUBLE nodes;

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
     * the interface (e.g. file reader) will provide external numbers of the
     * bounding nodes, so the constructor will adjust them to internal numbers.
     */
    VF_INT cells_vertices;

    //! Flag to indicate whether this is a submesh.
    bool submesh;

    //________________________________________//
    // End of private data of class TET_Builder.

 public:

    //! Explicit constructor, templated on Interface Type.
    template<class IT>
    explicit TET_Builder(SP<IT>);

    //! Build a TET_Mesh, using TET_Builder's private data.
    SP<TET_Mesh> build_Mesh();

};  // end class TET_Builder

//___________________________________________________________________________//
/*!
 * \brief           Non-inline, templated TET_Builder constructor.
 * \param interface Smart pointer to an instance of an Interface Type.
 */
template<class IT>
TET_Builder::TET_Builder(SP<IT> interface)
{
    Require (interface);

    // Get data arrays from the interface for TET_Mesh objects.
    coord_system   = interface->get_coord_system();
    nodes          = interface->get_nodes();
    parent         = interface->get_parent();
    cells_vertices = interface->get_cells_nodes();
    submesh        = interface->get_submesh();

    Ensure (coord_system == "xyz" || coord_system == "XYZ");

    // To be properly constructed, cells_vertices[][] must contain the
    // internal numbers of the four bounding vertices of the cells.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
        for (int ver_ = 0 ; ver_ < FOUR ; ver_++)
        {
            cells_vertices[cell_][ver_] -= 1;
            Ensure ( cells_vertices[cell_][ver_] >= 0 );
            Ensure ( cells_vertices[cell_][ver_] < nodes[0].size() );
        }

    // To be properly constructed, cells_vertices[][] must be in a predefined
    // order, so that later sense-of-surface operations will be predictable.

    for (int cell_ = 0 ; cell_ < cells_vertices.size() ; cell_++)
    {
        int v0 = cells_vertices[cell_][0];
        int v1 = cells_vertices[cell_][1];
        int v2 = cells_vertices[cell_][2];
        int v3 = cells_vertices[cell_][3];

        ThreeVector R0(nodes[0][v0],nodes[1][v0],nodes[2][v0]);
        ThreeVector R1(nodes[0][v1],nodes[1][v1],nodes[2][v1]);
        ThreeVector R2(nodes[0][v2],nodes[1][v2],nodes[2][v2]);
        ThreeVector R3(nodes[0][v3],nodes[1][v3],nodes[2][v3]);

        double disc = (R1 - R0).dot((R2 - R1).cross(R3 - R2));

        Check (disc != 0.0);

        if (disc < 0.0)
        {
            int temp = cells_vertices[cell_][2];
            cells_vertices[cell_][2] = cells_vertices[cell_][3];
            cells_vertices[cell_][3] = temp;
        }
    }

}   // end TET_Builder::TET_Builder(SP<IT>)

}   // end namespace rtt_mc

#endif  // __mc_TET_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/TET_Builder.hh
//---------------------------------------------------------------------------//
