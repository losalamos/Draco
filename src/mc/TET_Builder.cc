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
#include <utility>
#include <algorithm>
#include "TET_Builder.hh"

#ifdef __sun
#include <ctime>
#endif

namespace rtt_mc
{

using std::cerr;  // FOR DEBUGGING and TIMING.

using std::endl;
using std::set;
using std::make_pair;

//! Build a TET_Mesh, using TET_Builder's private data.
TET_Builder::SP_Mesh TET_Builder::build_Mesh()
{
    // Create first constructor argument: coord.
    TET_Builder::SP_Coord_sys coord(new XYZCoord_sys);

    // Create second constructor argument: layout.
    Layout layout(cells_vertices.size());
    for (int c = 1; c <= cells_vertices.size(); ++c)
        layout.set_size(c,FOUR);

    typedef std::map<SetInt, std::pair<int, int> > FaceMap;

    FaceMap faces;

#ifdef __sun
    clock_t t1 = clock();
#endif

    for (int cell = 0 ; cell < cells_vertices.size() ; ++cell)
    {
        FaceMap::iterator pos;

        SetInt face_0;
        face_0.insert(parent[cells_vertices[cell][1]]);
        face_0.insert(parent[cells_vertices[cell][2]]);
        face_0.insert(parent[cells_vertices[cell][3]]);

        pos = faces.find(face_0);
        if (pos == faces.end())
            faces.insert(make_pair(face_0, make_pair(cell+1, 1)));
        else
        {
            layout((pos->second).first, (pos->second).second) = cell +1;
            layout(cell+1, 1) = (pos->second).first;
        }

        SetInt face_1;
        face_1.insert(parent[cells_vertices[cell][2]]);
        face_1.insert(parent[cells_vertices[cell][3]]);
        face_1.insert(parent[cells_vertices[cell][0]]);

        pos = faces.find(face_1);
        if (pos == faces.end())
            faces.insert(make_pair(face_1, make_pair(cell+1, 2)));
        else
        {
            layout((pos->second).first, (pos->second).second) = cell +1;
            layout(cell+1, 2) = (pos->second).first;
        }

        SetInt face_2;
        face_2.insert(parent[cells_vertices[cell][3]]);
        face_2.insert(parent[cells_vertices[cell][0]]);
        face_2.insert(parent[cells_vertices[cell][1]]);

        pos = faces.find(face_2);
        if (pos == faces.end())
            faces.insert(make_pair(face_2, make_pair(cell+1, 3)));
        else
        {
            layout((pos->second).first, (pos->second).second) = cell +1;
            layout(cell+1, 3) = (pos->second).first;
        }

        SetInt face_3;
        face_3.insert(parent[cells_vertices[cell][0]]);
        face_3.insert(parent[cells_vertices[cell][1]]);
        face_3.insert(parent[cells_vertices[cell][2]]);

        pos = faces.find(face_3);
        if (pos == faces.end())
            faces.insert(make_pair(face_3, make_pair(cell+1, 4)));
        else
        {
            layout((pos->second).first, (pos->second).second) = cell +1;
            layout(cell+1, 4) = (pos->second).first;
        }
    }

#ifdef __sun
    clock_t t2 = clock();
    cerr << "Layout builder time = "
         << static_cast<double>(t2 - t1)/static_cast<double>(CLOCKS_PER_SEC)
         << endl;
#endif

    // Create third constructor argument: vertex_vector.
    sf_ThreeVector vertex_vector;
    for (int v_ = 0 ; v_ < node_coords.size() ; v_++)
        vertex_vector.push_back(ThreeVector(node_coords[v_][0],
                                node_coords[v_][1],node_coords[v_][2]));

    // Remaining constructor arguments (cells_vertices, submesh) already made.

    // Instantiate and return Smart Pointer to the new TET_Mesh.
    SP_Mesh mesh_ptr(new TET_Mesh(title, coord, layout,
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
