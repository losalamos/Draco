//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Builder.cc
 * \author H. Grady Hughes
 * \date   Mon Jan 31 13:23:43 MST 2000
 * \brief  Tetrahedral mesh builder class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <set>
#include <algorithm>
#include "TET_Builder.hh"

namespace rtt_mc
{

using std::endl;
using std::set;

//! Build a TET_Mesh, using TET_Builder's private data.
rtt_dsxx::SP<TET_Mesh> TET_Builder::build_Mesh()
{
    // Create first constructor argument: coord.
    rtt_dsxx::SP<Coord_sys> coord(new XYZCoord_sys);

    // Create second constructor argument: layout.
    Layout layout(cells_vertices.size());
    for (int c = 1; c <= cells_vertices.size(); ++c)
        layout.set_size(c,FOUR);

    int see_intersect[FOUR];         // Could actually be THREE.
    SetInt Face;

    for (int L_cell = 0; L_cell < cells_vertices.size() - 1; ++L_cell)
    {
        // Load L_total with all four vertices of L_cell.
        SetInt L_total;
        for (int v = 0; v < FOUR; ++v)
            L_total.insert(parent[cells_vertices[L_cell][v]]);

        for (int R_cell = L_cell + 1; R_cell < cells_vertices.size(); ++R_cell)
        {
            // Load R_total with all four vertices of R_cell.
            SetInt R_total;
            for (int v = 0; v < FOUR; ++v)
                R_total.insert(parent[cells_vertices[R_cell][v]]);

            int* end_intersect = std::set_intersection(L_total.begin(),
                    L_total.end(),R_total.begin(),R_total.end(),see_intersect);
            int N = end_intersect - see_intersect;
            Check (N >= 0 && N < FOUR);

            if (N == THREE)                         // Common face detected.
            {
                int face_found = -1;                // Work on L_cell neighbor.
                for (int v = 0; v < FOUR && face_found < 0; ++v)
                {
                    for (int f = 1 ; f <= THREE ; ++f)
                        Face.insert(parent[cells_vertices[L_cell][(f+v) % FOUR]]);

                    if (std::includes(R_total.begin(),R_total.end(),
                                      Face.begin(),Face.end()))
                        face_found = v;

                    Face.erase(Face.begin(),Face.end());
                }
                Check (face_found >= 0 && face_found < FOUR);

                layout(L_cell + 1, face_found + 1) = R_cell + 1;

                face_found = -1;                    // Work on R_cell neighbor.
                for (int v = 0; v < FOUR && face_found < 0; ++v)
                {
                    for (int f = 1 ; f <= THREE ; ++f)
                        Face.insert(parent[cells_vertices[R_cell][(f+v) % FOUR]]);

                    if (std::includes(L_total.begin(),L_total.end(),
                                      Face.begin(),Face.end()))
                        face_found = v;

                    Face.erase(Face.begin(),Face.end());
                }
                Check (face_found >= 0 && face_found < FOUR);

                layout(R_cell + 1, face_found + 1) = L_cell + 1;
            }
        }
    }

    // Create third constructor argument: vertex_vector.
    SF_THREEVECTOR vertex_vector;
    for (int v_ = 0 ; v_ < node_coords.size() ; v_++)
        vertex_vector.push_back(ThreeVector(node_coords[v_][0],
                                node_coords[v_][1],node_coords[v_][2]));

    // Remaining constructor arguments (cells_vertices, submesh) already made.

    // Instantiate and return Smart Pointer to the new TET_Mesh.
    rtt_dsxx::SP<TET_Mesh> mesh_ptr(new TET_Mesh(title, coord, layout,
        vertex_vector, node_coord_units, node_sets, side_sets, cell_sets,
        sides_vertices, cells_vertices, submesh));

    return mesh_ptr;

}   // end TET_Builder::build_Mesh()

//! Print the node_sets.
void TET_Builder::print_node_sets(std::ostream &output) const
{
    output << "NODE SETS\n" << endl;
    for (MAP_String_SetInt::const_iterator flag = node_sets.begin() ;
            flag != node_sets.end() ; flag++)
    {
        output << (*flag).first << "\n";
        int nnode = (*flag).second.size();
        output << "  " << nnode << " flagged node" <<
            (nnode == 1 ? "." : "s.") << "\n";
        if (nnode > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "    " << *i << "\n";
        output << endl;
    }
}   // end TET_Builder::print_node_sets(std::ostream &)

//! Print the element_sets.
void TET_Builder::print_element_sets(std::ostream &output) const
{
    output << "ELEMENT SETS\n" << endl;
    for (MAP_String_SetInt::const_iterator flag = element_sets.begin() ;
            flag != element_sets.end() ; flag++)
    {
        output << (*flag).first << "\n";
        int nelem = (*flag).second.size();
        output << "  " << nelem << " flagged element" <<
            (nelem == 1 ? "." : "s.") << "\n";
        if (nelem > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "    " << *i << "\n";
        output << endl;
    }
}   // end TET_Builder::print_element_sets(std::ostream &)

//! Print the side_sets.
void TET_Builder::print_side_sets(std::ostream &output) const
{
    output << "SIDE SETS\n" << endl;
    for (MAP_String_SetInt::const_iterator flag = side_sets.begin() ;
            flag != side_sets.end() ; flag++)
    {
        output << (*flag).first << "\n";
        int nside = (*flag).second.size();
        output << "  " << nside << " flagged side" <<
            (nside == 1 ? "." : "s.") << "\n";
        if (nside > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "    " << *i << "\n";
        output << endl;
    }
}   // end TET_Builder::print_side_sets(std::ostream &)

//! Print the cell_sets.
void TET_Builder::print_cell_sets(std::ostream &output) const
{
    output << "CELL SETS\n" << endl;
    for (MAP_String_SetInt::const_iterator flag = cell_sets.begin() ;
            flag != cell_sets.end() ; flag++)
    {
        output << (*flag).first << "\n";
        int ncell = (*flag).second.size();
        output << "  " << ncell << " flagged cell" <<
            (ncell == 1 ? "." : "s.") << "\n";
        if (ncell > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "    " << *i << "\n";
        output << endl;
    }
}   // end TET_Builder::print_cell_sets(std::ostream &)


}   // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of mc/TET_Builder.cc
//---------------------------------------------------------------------------//
