//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Builder.cc
 * \author H. Grady Hughes
 * \date   Mon Jan 31 13:23:43 MST 2000
 * \brief  Tetrahedral mesh builder class implementation file.
 */
//---------------------------------------------------------------------------//

#include <set>
#include <algorithm>
#include "TET_Builder.hh"

namespace rtt_mc
{

//! Build a TET_Mesh, using TET_Builder's private data.
rtt_dsxx::SP<TET_Mesh> TET_Builder::build_Mesh()
{
    // Create first constructor argument: coord.
    rtt_dsxx::SP<Coord_sys> coord(new XYZCoord_sys);

    // Create second constructor argument: layout.
    Layout layout(cells_vertices.size());
    for (int c = 1; c <= cells_vertices.size(); c++)
        layout.set_size(c,FOUR);

    std::vector< std::set<int> > L_face;
    L_face.resize(FOUR);

    std::vector< std::set<int> > R_face;
    R_face.resize(FOUR);

    std::set<int> R_total;

    for (int L_cell = 0; L_cell < cells_vertices.size() - 1; L_cell++)
    {
        // Load L_face[0] with vertices #1,#2,#3; L_face[1] with #2,#3,#0; etc.
        for (int f = 0; f < FOUR; f++)
            for (int v = 1; v <= THREE; v++)
                L_face[f].insert(parent[cells_vertices[L_cell][(f+v) % FOUR]]);

        for (int R_cell = L_cell + 1; R_cell < cells_vertices.size(); R_cell++)
        {
            // Load R_total with all four vertices of R_cell.
            for (int v = 0; v < FOUR; v++)
                R_total.insert(parent[cells_vertices[R_cell][v]]);

            // Set L_face_found if match is found.
            int L_face_found = -1;
            for (int f = 0; f < FOUR && L_face_found < 0; f++)
                if (std::includes(R_total.begin(),R_total.end(),
                             L_face[f].begin(),L_face[f].end()))
                    L_face_found = f;

            if (L_face_found >=0)
            {
                layout(L_cell + 1, L_face_found + 1) = R_cell + 1;

                // Load R_face[0] with vertices #1,#2,#3; etc.
                for (int f = 0; f < FOUR; f++)
                    for (int v = 1; v <= THREE; v++)
                        R_face[f].insert(parent[cells_vertices
                                          [R_cell][(f+v) % FOUR]]);

                // Set R_face_found when match is found.
                int R_face_found = -1;
                for (int f = 0; f < FOUR && R_face_found < 0; f++)
                    if (R_face[f] == L_face[L_face_found])
                        R_face_found = f;

                Check (R_face_found >= 0);
                layout(R_cell + 1, R_face_found + 1) = L_cell + 1;

                // Empty out R_face[0], R_face[1], R_face[2], and R_face[3].
                for (int f = 0; f < FOUR; f++)
                    R_face[f].erase(R_face[f].begin(),R_face[f].end());
            }

            // Empty out R_total.
            R_total.erase(R_total.begin(),R_total.end());
        }

        // Empty out L_face[0], L_face[1], L_face[2], and L_face[3].
        for (int f = 0; f < FOUR; f++)
            L_face[f].erase(L_face[f].begin(),L_face[f].end());
    }


    // Create third constructor argument: vertex_vector.
    SF_THREEVECTOR vertex_vector;
    for (int v_ = 0 ; v_ < nodes_coords.size() ; v_++)
        vertex_vector.push_back(ThreeVector(nodes_coords[v_][0],
                                nodes_coords[v_][1],nodes_coords[v_][2]));

    // Remaining constructor arguments (cells_vertices, submesh) already made.

    // Instantiate and return Smart Pointer to the new TET_Mesh.
    rtt_dsxx::SP<TET_Mesh> mesh_ptr(new
        TET_Mesh(coord, layout, vertex_vector, cells_vertices, submesh));

    return mesh_ptr;

}   // end TET_Builder::build_Mesh()

}   // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of mc/TET_Builder.cc
//---------------------------------------------------------------------------//
