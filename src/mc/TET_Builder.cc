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

    SF_INT neighbors_found(cells_vertices.size(),0);
    VF_INT cells_vert_parent(cells_vertices);
    for (int c = 0 ; c < cells_vertices.size() ; ++c)
        for (int v = 0 ; v < cells_vertices[c].size() ; ++v)
            cells_vert_parent[c][v] = parent[cells_vertices[c][v]];

    for (int L_cell = 0 ; L_cell < cells_vertices.size() - 1 ; ++L_cell)
    {
        if (neighbors_found[L_cell] == FOUR)
            continue;

        for (int R_cell = L_cell + 1; R_cell < cells_vertices.size(); ++R_cell)
        {
            if (neighbors_found[R_cell] == FOUR)
                continue;

            int count = 0;
            SF_INT L_hits(FOUR, -1);

            for (int lv = 0 ; lv < FOUR ; ++lv)
            {
                SF_INT::iterator R_itr = std::find(
                            cells_vert_parent[R_cell].begin(),
                            cells_vert_parent[R_cell].end(),
                            cells_vert_parent[L_cell][lv]);

                if (R_itr != cells_vert_parent[R_cell].end())
                {
                    ++count;
                    L_hits[lv] = R_itr - cells_vert_parent[R_cell].begin();
                }
            }
            if (count == THREE)
            {
                ++neighbors_found[L_cell];
                SF_INT::iterator L_itr = std::find(L_hits.begin(),
                                                    L_hits.end(), -1);
                int v = L_itr - L_hits.begin();
                layout(L_cell + 1, v + 1) = R_cell + 1;

                for (int rv = 0 ; rv < FOUR ; ++rv)
                    if (std::find(L_hits.begin(), L_hits.end(), rv)
                                               == L_hits.end())
                    {
                        ++neighbors_found[R_cell];
                        layout(R_cell + 1, rv + 1) = L_cell + 1;
                        break;
                    }
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
