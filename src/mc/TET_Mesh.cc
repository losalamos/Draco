//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Mesh.cc
 * \author H. Grady Hughes
 * \date   Tue Jan 18 10:15:43 MST 2000
 * \brief  Tetrahedral mesh class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TET_Mesh.hh"
#include "Constants.hh"
#include <iomanip>

namespace rtt_mc
{

//___________________________________________________________________________//
/*!
 * \brief                 TET_Mesh constructor.
 * \param coord_          Smart pointer to base class of derived XYZ system.
 * \param layout_         Layout of mesh.
 * \param vertex_vector_  ThreeVector vertices in a list of all vertices.
 * \param cells_vertices_ Internal identifiers of the four vertices of cells.
 * \param sides_vertices_ Internal identifiers of the three vertices of sides.
 * \param submesh_        Submesh indicator flag.
 */
TET_Mesh::TET_Mesh(std::string & title_, SP_Coord_sys coord_,
    Layout & layout_, sf_ThreeVector & vertex_vector_,
    std::string & node_coord_units_, MAP_String_SetInt & node_sets_,
    MAP_String_SetInt & side_sets_, MAP_String_SetInt & cell_sets_,
    vf_int & sides_vertices_, vf_int & cells_vertices_, bool submesh_)
    : title(title_), coord(coord_), layout(layout_),
      vertex_vector(vertex_vector_), node_coord_units(node_coord_units_),
      node_sets(node_sets_), side_sets(side_sets_), cell_sets(cell_sets_),
      sides_vertices(sides_vertices_), cells_vertices(cells_vertices_),
      submesh(submesh_)
{
    // For a TET_Mesh, there must be an XYZ coordinate system.
    Check (coord);
    Check (THREE == coord->get_dim());
    Check (std::string("xyz") == coord->get_Coord());

    // Confirm consistency of layout and cells_vertices.
    Check (layout.num_cells() == cells_vertices.size());
}

//___________________________________________________________________________//
/*!
 * \brief       Unnormalized outward cross-product of two edges of a face.
 * \param cell_ Internal number of cell.
 * \param face_ Internal identifier for the face (0, 1, 2, or 3).
 * \return      The unnormalized outward cross-product vector.
 *
 * Private function.
 *
 * The calling function will have checked the validity of the external
 * values of "cell" and "face" and will have adjusted to cell_ and face_
 * before calling this private function.
 *
 * The magnitude of this vector is twice the face surface area.
 *
 * This function is used by get_normal().
 */
const ThreeVector TET_Mesh::get_outward_cross(int cell_, int face_) const
{
    int head1, tail1, head2, tail2;

    if (face_ == 0)
    {
        head1 = tail2 = cells_vertices[cell_][2];
        tail1 = cells_vertices[cell_][1];
        head2 = cells_vertices[cell_][3];
    }
    else if (face_ == 1)
    {
        head1 = cells_vertices[cell_][2];
        tail1 = tail2 = cells_vertices[cell_][3];
        head2 = cells_vertices[cell_][0];
    }
    else if (face_ == 2)
    {
        head1 = tail2 = cells_vertices[cell_][0];
        tail1 = cells_vertices[cell_][3];
        head2 = cells_vertices[cell_][1];
    }
    else
    {
        head1 = cells_vertices[cell_][0];
        tail1 = tail2 = cells_vertices[cell_][1];
        head2 = cells_vertices[cell_][2];
    }

    return (vertex_vector[head1] - vertex_vector[tail1]).cross(
                vertex_vector[head2] - vertex_vector[tail2]);


}   // end TET_Mesh::get_outward_cross(int,int)

//___________________________________________________________________________//
/*!
 * \brief       Unnormalized inward cross-product of two edges of a face.
 * \param cell_ Internal number of cell.
 * \param face_ Internal identifier for the face (0, 1, 2, or 3).
 * \return      The unnormalized inward cross-product vector.
 *
 * Private function.
 *
 * The calling function will have checked the validity of the external
 * values of "cell" and "face" and will have adjusted to cell_ and face_
 * before calling this private function.
 *
 * The magnitude of this vector is twice the face surface area.
 *
 * This function is used by get_normal_in(), by volume(), and by face_area().
 */
const ThreeVector TET_Mesh::get_inward_cross(int cell_, int face_) const
{
    int head1, tail1, head2, tail2;

    if (face_ == 0)
    {
        head1 = cells_vertices[cell_][1];
        tail1 = tail2 = cells_vertices[cell_][2];
        head2 = cells_vertices[cell_][3];
    }
    else if (face_ == 1)
    {
        head1 = tail2 = cells_vertices[cell_][3];
        tail1 = cells_vertices[cell_][2];
        head2 = cells_vertices[cell_][0];
    }
    else if (face_ == 2)
    {
        head1 = cells_vertices[cell_][3];
        tail1 = tail2 = cells_vertices[cell_][0];
        head2 = cells_vertices[cell_][1];
    }
    else
    {
        head1 = tail2 = cells_vertices[cell_][1];
        tail1 = cells_vertices[cell_][0];
        head2 = cells_vertices[cell_][2];
    }

    return (vertex_vector[head1] - vertex_vector[tail1]).cross(
                vertex_vector[head2] - vertex_vector[tail2]);

}   // end TET_Mesh::get_inward_cross(int,int)

//___________________________________________________________________________//
/*!
 * \brief          Calculate the barycentric coordinates of a point in a cell.
 * \param position XYZ-position of the point as STL vector.
 * \param cell_    Internal number of cell.
 * \return         The normalized barycentric coordinates of "position".
 *
 * Private function.
 *
 * The calling function will have checked the validity of the external value of
 * "cell" and will have adjusted to cell_ before calling this private function.
 *
 * If the position is not actually within cell_, a Check() will fail.
 *
 * The barycentric coordinates will be normalized and ordered as the vertices
 * (and corresponding faces) are ordered - 0, 1, 2, 3.
 */
const TET_Mesh::sf_double TET_Mesh::get_barycentric_coords(
        const sf_double &position, int cell_) const
{
    sf_double b_coords(FOUR);
    ThreeVector R(position);
    double summ = 0.;

    for (int f_ = 0 ; f_ < FOUR ; f_++)
    {
        ThreeVector N = get_inward_cross(cell_, f_);
        int v_ = (f_ + 1) % FOUR;  // any vertex on the face.
        int base = cells_vertices[cell_][v_];

        b_coords[f_] = N.dot(R - vertex_vector[base]);
        Check (b_coords[f_] > 0.);
        summ += b_coords[f_];
    }

    for (int f_ = 0 ; f_ < FOUR ; f_++)
        b_coords[f_] /= summ;

    return b_coords;

}   // end TET_Mesh::get_barycentric_coords(const sf_double &,int)

//___________________________________________________________________________//
/*!
 * \brief          Determine whether a given position is inside a given cell.
 * \param position XYZ-position as STL vector.
 * \param cell     External number of cell.
 * \return         True if position is strictly within the cell, else false.
 *
 * The external number "cell" is checked for validity, then converted to
 * cell_ = cell - 1 for internal use.
 */
bool TET_Mesh::in_open_cell(const sf_double &position, int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);

    for (int f_ = 0 ; f_ < FOUR ; f_++)
    {
        ThreeVector N = get_outward_cross(cell_, f_);
        int v_ = (f_ + 1) % FOUR;  // any vertex on the face.

        if ( N.dot(vertex_vector[cells_vertices[cell_][v_]] - XYZ) <= 0.0 )
            return false;
    }

    // If position is "inside" every face, it is inside the cell.
    return true;

}   // end TET_Mesh::in_open_cell(const sf_double &, int)

//___________________________________________________________________________//
/*!
 * \brief          Determine if a position is in or on the boundary of a cell.
 * \param position XYZ-position as STL vector.
 * \param cell     External number of cell.
 * \return         True if position is in the cell or on the boundary.
 *
 * The external number "cell" is checked for validity, then converted to
 * cell_ = cell - 1 for internal use.
 */
bool TET_Mesh::in_closed_cell(const sf_double &position, int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);

    for (int f_ = 0 ; f_ < FOUR ; f_++)
    {
        ThreeVector N = get_outward_cross(cell_, f_);
        int v_ = (f_ + 1) % FOUR;  // any vertex on the face.

        if ( N.dot(vertex_vector[cells_vertices[cell_][v_]] - XYZ) < 0.0 )
            return false;
    }
    // If position is "inside-or-on" every face, it is in the closed cell.
    return true;

}   // end TET_Mesh::in_closed_cell(const sf_double &, int)

//___________________________________________________________________________//
/*!
 * \brief          Find the distance to the boundary and the encountered face.
 * \param position XYZ-position of particle.
 * \param omega    XYZ-direction of particle.
 * \param cell     External number of current cell (containing position).
 * \param face     (returned reference) External number of encountered face.
 * \return         The distance to the encountered boundary face.
 *
 * The external number "cell" is checked for validity, then converted to
 * cell_ = cell - 1 for internal use.
 *
 * The internal identifier face_ will be discovered, and converted to
 * face = face_ + 1 for return, so that face == 1, 2, 3, or 4.
 *
 * Future considerations:
 *
 * 1. That "position" is actually in "cell" is not checked.
 *
 * 2. The normalization of omega is not checked.
 *
 * 3. The function does not check for encounters with cell edges or vertices.
 *
 * 4. No consideration of machine round-off error is made.
 *
 * 5. For notation, both position and omega are converted to ThreeVectors.
 * This could be hand-coded for efficiency if necessary.
 */
double TET_Mesh::get_db(const sf_double &position, const sf_double &omega,
                        int cell, int &face) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);
    ThreeVector UVW(omega);
    sf_double dist(FOUR);

    for (int f_ = 0 ; f_ < FOUR ; f_++)
    {
        ThreeVector N = get_outward_cross(cell_, f_);
        double denom = N.dot(UVW);
        int v = (f_ + 1) % FOUR;  // any vertex on the face.

        if (denom > 0.0)
            dist[f_] = N.dot(vertex_vector[cells_vertices[cell_][v]] - XYZ)
                       /denom;
        else
            dist[f_] = global::huge;
    }

    sf_double::iterator itor = std::min_element(dist.begin(),dist.end());

    double dist_min = *itor;
    Ensure ( dist_min > 0.0 && dist_min < global::huge );

    int face_ = itor - dist.begin();
    face = face_ + 1;
    Ensure ( face >= 1 && face <= FOUR );

    return *itor;

}   // end TET_Mesh::get_db(const sf_double&,const sf_double&,int,int&)

//___________________________________________________________________________//
/*!
 * \brief      Calculate the outward normal vector for a given cell and face.
 * \param cell External number of cell.
 * \param face External identifier for the face (1, 2, 3, or 4).
 * \return     The normalized (unit-magnitude) outward normal vector.
 *
 * The external numbers "cell" and "face" are checked for validity, and the
 * internal numbers "cell-1" and "face-1" are passed to get_outward_cross().
 */
const TET_Mesh::sf_double TET_Mesh::get_normal(int cell, int face) const
{
    Valid(cell, face);
    ThreeVector N = get_outward_cross(cell-1, face-1);

    N.normalize();

    return N.convert();

}   // end TET_Mesh::get_normal(int,int)

//___________________________________________________________________________//
/*!
 * \brief      Calculate the inward normal vector for a given cell and face.
 * \param cell External number of cell.
 * \param face External identifier for the face (1, 2, 3, or 4).
 * \return     The normalized (unit-magnitude) inward normal vector.
 *
 * The external numbers "cell" and "face" are checked for validity, and the
 * internal numbers "cell-1" and "face-1" are passed to get_outward_cross().
 */
const TET_Mesh::sf_double TET_Mesh::get_normal_in(int cell, int face) const
{
    Valid(cell, face);
    ThreeVector N = get_inward_cross(cell-1, face-1);

    N.normalize();

    return N.convert();

}   // end TET_Mesh::get_normal_in(int,int)

//___________________________________________________________________________//
/*!
 * \brief      Calculate the volume of a cell.
 * \param cell External number of cell.
 * \return     The volume of the cell.
 *
 * The external number "cell" is checked for validity.  Then the internal
 * number cell_ = cell - 1 is passed to get_inward_cross() and used locally.
 */
double TET_Mesh::volume(int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;
    ThreeVector N = get_inward_cross(cell_, 0);

    int v0 = cells_vertices[cell_][0];
    int v1 = cells_vertices[cell_][1];

    double vol = N.dot(vertex_vector[v0] - vertex_vector[v1])/6.0;

    Ensure ( vol > 0.0 );

    return vol;

}   // end TET_Mesh::volume(int)

//___________________________________________________________________________//
/*!
 * \brief      Calculate the surface area of a given face of a given cell.
 * \param cell External number of cell.
 * \param face External identifier for the face (1, 2, 3, or 4).
 * \return     The surface area.
 *
 * The external numbers "cell" and "face" are checked for validity, and the
 * internal numbers "cell-1" and "face-1" are passed to get_inward_cross().
 */
double TET_Mesh::face_area(int cell, int face) const
{
    Valid(cell, face);
    double area = 0.5*get_inward_cross(cell-1, face-1).get_norm();

    Ensure ( area > 0.0 );

    return area;

}   // end TET_Mesh::face_area(int,int)

//___________________________________________________________________________//
/*!
 * \brief      Get coordinates of the four vertices bounding a given cell.
 * \param cell External number of cell.
 * \return     Vector_field[dim#][vertex#] == coordinate along dim#-axis.
 *
 * dim# = 0, 1, 2 <----> for axis X, Y, Z.
 *
 * vertex# = 0, 1, 2, 3.
 */
const TET_Mesh::vf_double TET_Mesh::get_vertices(int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;
    // Could check that (cells_vertices[cell_].size() == FOUR);

    TET_Mesh::vf_double ret_vert(THREE);

    for (int cv = 0 ; cv < FOUR ; cv++)
    {
        int v = cells_vertices[cell_][cv];
        ret_vert[0].push_back(vertex_vector[v].get_x());
        ret_vert[1].push_back(vertex_vector[v].get_y());
        ret_vert[2].push_back(vertex_vector[v].get_z());
    }

    return ret_vert;

}   // end TET_Mesh::get_vertices(int)

//___________________________________________________________________________//
/*!
 * \brief      Get coordinates of the three vertices bounding a given face.
 * \param cell External number of cell.
 * \param face External identifier for the face (1, 2, 3, or 4).
 * \return     Vector_field[dim#][vertex#] == coordinate along dim#-axis.
 *
 * dim# = 0, 1, 2 <----> for axis X, Y, Z.
 *
 * vertex# = 0, 1, 2, but representing three out of four of (0, 1, 2, 3),
 * since each face has 3 vertices.  For no particular reason, these are
 * chosen so that right-handed circulation in the order specified would
 * result in an outward direction.
 */
const TET_Mesh::vf_double TET_Mesh::get_vertices(int cell, int face) const
{
    Valid(cell, face);
    int cell_ = cell - 1;
    int face_ = face - 1;

    vf_double ret_vert(THREE);

    int fv[THREE];

    if (face_ == 0)
    {
        fv[0] = 1; fv[1] = 2; fv[2] = 3;
    }
    else if (face_ == 1)
    {
        fv[0] = 2; fv[1] = 0; fv[2] = 3;
    }
    else if (face_ == 2)
    {
        fv[0] = 3; fv[1] = 0; fv[2] = 1;
    }
    else
    {
        fv[0] = 0; fv[1] = 2; fv[2] = 1;
    }

    for (int i = 0 ; i < THREE ; i++)
    {
        int v = cells_vertices[cell_][fv[i]];
        ret_vert[0].push_back(vertex_vector[v].get_x());
        ret_vert[1].push_back(vertex_vector[v].get_y());
        ret_vert[2].push_back(vertex_vector[v].get_z());
    }

    return ret_vert;

}   // end TET_Mesh::get_vertices(int,int)

//___________________________________________________________________________//
/*!
 * \brief        Sample position uniformly within a given cell.
 * \param cell   External number of cell.
 * \param random The random-number generator, used as random.ran()
 * \return       Scalar_field[dim#] == sampled coordinate along dim#-axis.
 *
 * This version samples uniformly for a position on the triangular base first,
 * then samples along the altitude, with the appropriate geometric PDF, for the
 * position in the cell.
 */
const TET_Mesh::sf_double TET_Mesh::sample_pos(int cell,
    rtt_rng::Sprng &random) const
{
    Valid(cell);
    int cell_ = cell - 1;

    int v0 = cells_vertices[cell_][0];
    int v1 = cells_vertices[cell_][1];
    int v2 = cells_vertices[cell_][2];
    int v3 = cells_vertices[cell_][3];

    ThreeVector B = sample_in_triangle(vertex_vector[v0],
                                       vertex_vector[v1],
                                       vertex_vector[v2],random);

    double fraction = std::pow(random.ran(),1.0/3.0);
    Check ( fraction > 0.0 && fraction < 1.0 );

    return rtt_mc::lin_comb(vertex_vector[v3],B,fraction).convert();

}   // end TET_Mesh::sample_pos(int,rtt_rng::Sprng &)

//___________________________________________________________________________//
/*!
 * \brief         Sample position in a given cell with a T**4 distribution.
 * \param cell    External number of cell.
 * \param random  The random-number generator, used as random.ran()
 * \param T4      Values of T**4, in the same order as the cell's vertices.
 * \return        Scalar_field[dim#] == sampled coordinate along dim#-axis.
 *
 * This function assumes that the temperatures are known on the four vertices
 * of the cell, and are adequately represented by linear interpolation in T**4.
 */
const TET_Mesh::sf_double TET_Mesh::sample_pos(int cell,
    rtt_rng::Sprng &random, const sf_double &T4) const
{
    Require (T4.size() == FOUR);
    Valid(cell);
    int cell_ = cell - 1;
    for (int i = 0 ; i < FOUR ; i++)
        Check (T4[i] >= 0.);

    sf_double::const_iterator im = std::max_element(T4.begin(),T4.end());
    double Tmax = *im;
    Check (Tmax >= 0.);
    double Tr;

    sf_double R(THREE);
    do
    {
        R = sample_pos(cell, random);  // Trial uniform sampling.
        sf_double B = get_barycentric_coords(R, cell_);
        Tr = B[0]*T4[0] + B[1]*T4[1] + B[2]*T4[2] + B[3]*T4[3];
    } while (Tr <= Tmax*random.ran());

    return R;

}   // end TET_Mesh::sample_pos(int,rtt_rng::Sprng &,const sf_double &)

//___________________________________________________________________________//
/*!
 * \brief        Sample position uniformly on a given face of a given cell.
 * \param cell   External number of cell.
 * \param face   External identifier for the face (1, 2, 3, or 4).
 * \param random The random-number generator, used as random.ran()
 * \return       Scalar_field[dim#] == sampled coordinate along dim#-axis.
 */
const TET_Mesh::sf_double TET_Mesh::sample_pos_on_face(int cell, int face,
    rtt_rng::Sprng &random) const
{
    Valid(cell, face);
    int cell_ = cell - 1;
    int face_ = face - 1;
    int va, vb, vc;

    // Identify the vertices belonging to the face. For no particular reason,
    // choose these so that right-handed circulation in the order specified
    // would result in an outward direction.
    if (face_ == 0)
    {
        va = 1; vb = 2; vc = 3;
    }
    else if (face_ == 1)
    {
        va = 2; vb = 0; vc = 3;
    }
    else if (face_ == 2)
    {
        va = 3; vb = 0; vc = 1;
    }
    else
    {
        va = 0; vb = 2; vc = 1;
    }

    return rtt_mc::sample_in_triangle(vertex_vector[cells_vertices[cell_][va]],
                              vertex_vector[cells_vertices[cell_][vb]],
                              vertex_vector[cells_vertices[cell_][vc]],
                              random).convert();

}   // end TET_Mesh::sample_pos_on_face(int,int,rtt_rng::Sprng &)

//___________________________________________________________________________//
/*!
 * \brief Overloaded operator== for design-by-contract.
 */
bool TET_Mesh::operator==(const TET_Mesh &rhs) const
{
// QUESTION: Should we check submesh ?

    // Verify the mesh titles.
    if (title != rhs.title)
        return false;

    // Verify that we have the same coordinate systems.
    if (coord->get_Coord() != rhs.coord->get_Coord())
        return false;

    // Verify that the Layouts are equal.
    if (layout != rhs.layout)
        return false;

    // Verify the coordinates of the vertices.
    if (vertex_vector != rhs.vertex_vector)
        return false;

    // Verify the coordinate units.
    if (node_coord_units != rhs.node_coord_units)
        return false;

    // Verify the node_sets.
    if (node_sets != rhs.node_sets)
        return false;

    // Verify the side_sets.
    if (side_sets != rhs.side_sets)
        return false;

    // Verify the cell_sets.
    if (cell_sets != rhs.cell_sets)
        return false;

    // Verify the identities of the vertices of each side.
    if (sides_vertices != rhs.sides_vertices)
        return false;

    // Verify the identities of the vertices of each cell.
    if (cells_vertices != rhs.cells_vertices)
        return false;

    // if we haven't returned, then the two meshes must be equal
    return true;
}   // end TET_Mesh::operator==(const TET_Mesh &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of mesh title.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_title(std::ostream &output) const
{
    for (int i = 0 ; i < title.size() ; i++)
        output << "_";
    output << "\n" << title << std::endl;
}   // end TET_Mesh::print_title(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of layout.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_layout(std::ostream &output) const
{
    output << "___ LAYOUT ___\n";
    output << "Cell:  Face ==> Neighbor;  Face ==> Neighbor; ...";
    for (int c = 0 ; c < layout.num_cells() ; c++)
    {
        output << "\n    cell " << c << ":";
        for (int f = 0 ; f < layout.num_faces(c+1) ; f++)
        {
            output << "  " << f << " ==> ";
            int n = layout(c+1,f+1);
            if (n == 0)
                output << "ext;";
            else
                output << n-1 << ";";
        }
    }
    output << std::endl;
}   // end TET_Mesh::print_layout(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of vertex_vector.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_vertex_vector(std::ostream &output) const
{
    output << "___ VERTEX VECTOR ___";
    std::ios_base::fmtflags prev_state = output.flags();
    std::streamsize prev_prec = output.precision();
    output << std::scientific << std::setprecision(5);

    for (int v = 0 ; v < vertex_vector.size() ; v++)
        output << "\n   " << v << ":  " << vertex_vector[v].get_x() << ", " <<
            vertex_vector[v].get_y() << ", " << vertex_vector[v].get_z();
    output << std::endl;

    output.flags(prev_state);
    output << std::setprecision(prev_prec);
}   // end TET_Mesh::print_vertex_vector(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of node_coord_units.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_node_coord_units(std::ostream &output) const
{
    output << "___ Node coordinate units are " << node_coord_units <<
              " ___" << std::endl;
}   // end TET_Mesh::print_node_coord_units(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of node_sets.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_node_sets(std::ostream &output) const
{
    output << "___ NODE SETS ___";
    for (MAP_String_SetInt::const_iterator flag = node_sets.begin() ;
            flag != node_sets.end() ; flag++)
    {
        output << "\n" << (*flag).first;
        int nnode = (*flag).second.size();
        output << "\n   " << nnode << " flagged node" <<
            (nnode == 1 ? "." : "s.");
        if (nnode > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "\n      " << *i;
    }
    output << std::endl;
}   // end TET_Mesh::print_node_sets(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of side_sets.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_side_sets(std::ostream &output) const
{
    output << "___ SIDE SETS ___";
    for (MAP_String_SetInt::const_iterator flag = side_sets.begin() ;
            flag != side_sets.end() ; flag++)
    {
        output << "\n" << (*flag).first;
        int nside = (*flag).second.size();
        output << "\n   " << nside << " flagged side" <<
            (nside == 1 ? "." : "s.");
        if (nside > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "\n      " << *i;
    }
    output << std::endl;
}   // end TET_Mesh::print_side_sets(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of cell_sets.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_cell_sets(std::ostream &output) const
{
    output << "___ CELL SETS ___";
    for (MAP_String_SetInt::const_iterator flag = cell_sets.begin() ;
            flag != cell_sets.end() ; flag++)
    {
        output << "\n" << (*flag).first;
        int ncell = (*flag).second.size();
        output << "\n   " << ncell << " flagged cell" <<
            (ncell == 1 ? "." : "s.");
        if (ncell > 0)
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                                   i != (*flag).second.end(); i++)
                output << "\n      " << *i;
    }
    output << std::endl;
}   // end TET_Mesh::print_cell_sets(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of sides_vertices.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_sides_vertices(std::ostream &output) const
{
    output << "___ SIDES VERTICES ___";
    for (int s = 0 ; s < sides_vertices.size() ; s++)
    {
        output << "\n    side " << s << ":";
        for (int v = 0 ; v < sides_vertices[s].size() ; v++)
            output << "   " << sides_vertices[s][v];
    }
    output << std::endl;
}   // end TET_Mesh::print_sides_vertices(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of cells_vertices.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_cells_vertices(std::ostream &output) const
{
    output << "___ CELLS VERTICES ___";
    for (int c = 0 ; c < cells_vertices.size() ; c++)
    {
        output << "\n    cell " << c << ":";
        for (int v = 0 ; v < cells_vertices[c].size() ; v++)
            output << "   " << cells_vertices[c][v];
    }
    output << std::endl;
}   // end TET_Mesh::print_cells_vertices(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Allow output of submesh status.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_submesh(std::ostream &output) const
{
    if (submesh)
        output << "___ This is a submesh. ___\n" <<
                  "--------------------------" << std::endl;
    else
        output << "___ This is a full mesh. ___\n" <<
                  "----------------------------" << std::endl;
}   // end TET_Mesh::print_submesh(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief        Output everything about the mesh.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::print_mesh(std::ostream &output) const
{
    print_title(output);
    print_layout(output);
    print_vertex_vector(output);
    print_node_coord_units(output);
    print_node_sets(output);
    print_side_sets(output);
    print_cell_sets(output);
    print_sides_vertices(output);
    print_cells_vertices(output);
    print_submesh(output);
}   // end TET_Mesh::print_mesh(std::ostream &)

//___________________________________________________________________________//
/*!
 * \brief          Find the cell comtaining a given position.
 * \param position XYZ-position as STL vector.
 * \return         External number of cell containing position.
 *
 * This function should only be called for full meshes, and for positions
 * that are strictly within some cell.  Positions outside the mesh or lying
 * on boundaries will violate an Ensure() test.
 *
 * For TET_Meshes, no binary search is available, so get_cell() will be slow.
 */
int TET_Mesh::get_cell(const sf_double &position) const
{
    Require (!submesh);

    int ret_cell = 0;

    for (int cell = 1 ; cell <= num_cells() ; cell++)
        if ( in_open_cell(position, cell) )
            ret_cell = cell;

    Ensure (ret_cell > 0);
    return ret_cell;

}   // end TET_Mesh::get_cell(const sf_double &)

//___________________________________________________________________________//
/*!
 * \brief          Find the minimum distance to boundary in the given cell.
 * \param position XYZ-position of particle, as STL vector.
 * \param cell     External number of current cell (containing position).
 * \return         The distance to the encountered boundary face.
 *
 * The external number "cell" is checked for validity, then converted to
 * cell_ = cell - 1 for internal use.  However, "cell" and the for-loop
 * index "face" are passed unmodified to get_normal(cell, face).
 *
 * It is not checked that "position" is actually in "cell".
 *
 */
double TET_Mesh::get_min_db(const sf_double &position, int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);
    sf_double dist(FOUR);

    for (int face = 1 ; face <= FOUR ; face++)
    {
        ThreeVector N(get_normal(cell, face));
        int v_ = face % FOUR;  // any vertex on the face.

        dist[face-1] =
            N.dot(vertex_vector[cells_vertices[cell_][v_]] - XYZ);
    }

    sf_double::iterator itor = std::min_element(dist.begin(),dist.end());

    double dist_min = *itor;
    Ensure ( dist_min > 0.0 && dist_min < global::huge );
    return dist_min;

}   // end TET_Mesh::get_min_db(const sf_double &, int)

//===========================================================================//
// Beginning of TET_Mesh::Pack member functions
//===========================================================================//

//___________________________________________________________________________//
/*!
 * \brief    TET_Mesh::Pack constructor.
 * \param ds Size of data array of doubles.
 * \param dd Pointer to data array of doubles.
 * \param is Size of data array of integers.
 * \param id Pointer to data array of integers.
 * \param cs Size of data array of characters.
 * \param cd Pointer to data array of characters.
 *
 * Note that after construction, the Pack object owns the pointed-to data,
 * and assumes the responsibility of deletion.
 */
TET_Mesh::Pack::Pack(int ds, double *dd, int is, int *id, int cs, char *cd)
    : dsize(ds), ddata(dd), isize(is), idata(id), csize(cs), cdata(cd)
{
}

//---------------------------------------------------------------------------//
//! TET_Mesh::Pack destructor.
TET_Mesh::Pack::~Pack()
{
    delete [] ddata;
    delete [] idata;
    delete [] cdata;
}

//---------------------------------------------------------------------------//
//! TET_Mesh::Pack copy constuctor

TET_Mesh::Pack::Pack(const Pack &rhs)
    : dsize(rhs.dsize), ddata(new double[dsize]),
      isize(rhs.isize), idata(new int[isize]),
      csize(rhs.csize), cdata(new char[csize])
{
    // copy double data
    for (int i = 0; i < dsize; ++i)
        ddata[i] = rhs.ddata[i];

    // copy int data
    for (int i = 0; i < isize; ++i)
        idata[i] = rhs.idata[i];

    // copy char data
    for (int i = 0; i < csize; ++i)
        cdata[i] = rhs.cdata[i];
}

//---------------------------------------------------------------------------//
//! Unpacking routine: member function of struct Pack.
TET_Mesh::SP_Mesh TET_Mesh::Pack::unpack() const
{
    int LenTitle = idata[0];          // length of title string
    int LenUnits = idata[1];          // length of coord_units string
    int NumVerts = idata[2];          // number of vertices
    int NumSides = idata[3];          // number of sides
    int NumCells = idata[4];          // number of cells
    int SubmeshFlag = idata[5];       // 1 for submesh, else 0
    int NumNodeSets = idata[6];       // number of node sets
    int NumSideSets = idata[7];       // number of side sets
    int NumCellSets = idata[8];       // number of cell sets

    // Counters for assigning data.
    int d_ctr  = 0;           // for doubles.
    int i_ctr  = 9;           // for integers (already saw 0..8).
    int c_ctr  = 0;           // for characters.

   // TET meshes are always XYZ.
    TET_Mesh::SP_Coord_sys coord_(new XYZCoord_sys());

    bool submesh_ = (SubmeshFlag == 1);

    std::string title_;
    for (int t = 0 ; t < LenTitle ; ++t)
        title_ += cdata[c_ctr++];

    std::string node_coord_units_;
    for (int u = 0 ; u < LenUnits ; ++u)
        node_coord_units_ += cdata[c_ctr++];

    std::vector<ThreeVector> vertex_vector_;
    for (int v = 0 ; v < NumVerts ; ++v)
    {
        double x = ddata[d_ctr++];
        double y = ddata[d_ctr++];
        double z = ddata[d_ctr++];
        vertex_vector_.push_back(ThreeVector(x,y,z));
    }

    Layout layout_;
    layout_.set_size(NumCells);
    for (int c = 0 ; c < NumCells ; ++c)
        layout_.set_size(c+1,FOUR);

    for (int c = 0 ; c < NumCells ; ++c)
        for (int f = 0 ; f < FOUR ; ++f)
            layout_(c+1,f+1) = idata[i_ctr++];

    vf_int cells_vertices_(NumCells);
    for (int c = 0 ; c < NumCells ; ++c)
        for (int v = 0 ; v < FOUR ; ++v)
            cells_vertices_[c].push_back(idata[i_ctr++]);

    vf_int sides_vertices_(NumSides);
    for (int s = 0 ; s < NumSides ; ++s)
        for (int v = 0 ; v < THREE ; ++v)
            sides_vertices_[s].push_back(idata[i_ctr++]);

    MAP_String_SetInt node_sets_;

    if (NumNodeSets > 0)
    {
        for (int n = 0 ; n < NumNodeSets ; ++n)
        {
            int slen = idata[i_ctr++];
            int mems = idata[i_ctr++];
            std::string set_id;

            for (int c = 0 ; c < slen ; ++c)
                set_id += cdata[c_ctr++];

            SetInt set_mems;
            for (int m = 0 ; m < mems ; ++m)
                set_mems.insert(idata[i_ctr++]);

            node_sets_[set_id] = set_mems;
        }
    } // End if (NumNodeSets > 0)

    MAP_String_SetInt side_sets_;

    if (NumSideSets > 0)
    {
        for (int s = 0 ; s < NumSideSets ; ++s)
        {
            int slen = idata[i_ctr++];
            int mems = idata[i_ctr++];
            std::string set_id;

            for (int c = 0 ; c < slen ; ++c)
                set_id += cdata[c_ctr++];

            SetInt set_mems;
            for (int m = 0 ; m < mems ; ++m)
                set_mems.insert(idata[i_ctr++]);

            side_sets_[set_id] = set_mems;
        }
    } // End if (NumSideSets > 0)

    MAP_String_SetInt cell_sets_;

    if (NumCellSets > 0)
    {
        for (int c = 0 ; c < NumCellSets ; ++c)
        {
            int slen = idata[i_ctr++];
            int mems = idata[i_ctr++];
            std::string set_id;

            for (int c = 0 ; c < slen ; ++c)
                set_id += cdata[c_ctr++];

            SetInt set_mems;
            for (int m = 0 ; m < mems ; ++m)
                set_mems.insert(idata[i_ctr++]);

            cell_sets_[set_id] = set_mems;
        }
    } // End if (NumCellSets > 0)

    SP_Mesh mesh_ptr(new TET_Mesh(title_, coord_, layout_,
        vertex_vector_, node_coord_units_, node_sets_, side_sets_, cell_sets_,
        sides_vertices_, cells_vertices_, submesh_));

    return mesh_ptr;
}   // end TET_Mesh::Pack::unpack()

/*!
 * \brief        Allow output of packed meshes.
 * \param output The standard output stream to recieve the report.
 */
void TET_Mesh::Pack::print_pack(std::ostream &output) const
{
    output << "Begin print of packed mesh.\n";
    output << dsize << " double(s).\n";
    for (int d = 0 ; d < dsize ; ++d)
        output << ddata[d] << "\n";

    output << isize << " integer(s).\n";
    for (int i = 0 ; i < isize ; ++i)
        output << idata[i] << "\n";

    output << csize << " character(s).\n";
    for (int c = 0 ; c < csize ; ++c)
        output << cdata[c] << "\n";

    output << "End print of packed mesh.\n";
}   // end TET_Mesh::Pack::print_pack(std::ostream &)

//===========================================================================//
// End of TET_Mesh::Pack member functions
//===========================================================================//

//___________________________________________________________________________//
/*!
 * \brief TET_Mesh member function to return a Pack object holding a TET_Mesh.
 *
 * This initial version uses STL vectors to build the three arrays, allocating
 * conventional arrays later.  This makes the counting easier, but is wasteful
 * of memory.  A later version will count first, then allocate space only once.
 */
TET_Mesh::Pack TET_Mesh::pack() const
{
    int LenTitle = title.length();             // length of title string
    int LenUnits = node_coord_units.length();  // length of coord_units string
    int NumVerts = vertex_vector.size();       // number of vertices
    int NumSides = sides_vertices.size();      // number of sides
    int NumCells = cells_vertices.size();      // number of cells
    int SubmeshFlag = (submesh ? 1 : 0);       // 1 for submesh, else 0
    int NumNodeSets = node_sets.size();        // number of node sets
    int NumSideSets = side_sets.size();        // number of side sets
    int NumCellSets = cell_sets.size();        // number of cell sets

    // Counters for assigning data.
    int d_ctr  = 0;           // for doubles.
    int i_ctr  = 0;           // for integers.
    int c_ctr  = 0;           // for characters.

    // Temporary containers for data.
    std::vector<double> d_vec;       // for doubles.
    std::vector<int> i_vec;          // for integers.
    std::vector<char> c_vec;         // for characters.

    i_vec.push_back(LenTitle); ++i_ctr;
    i_vec.push_back(LenUnits); ++i_ctr;
    i_vec.push_back(NumVerts); ++i_ctr;
    i_vec.push_back(NumSides); ++i_ctr;
    i_vec.push_back(NumCells); ++i_ctr;
    i_vec.push_back(SubmeshFlag); ++i_ctr;
    i_vec.push_back(NumNodeSets); ++i_ctr;
    i_vec.push_back(NumSideSets); ++i_ctr;
    i_vec.push_back(NumCellSets); ++i_ctr;

    for (int c = 0 ; c < LenTitle ; ++c)
    {
        c_vec.push_back(title[c]); ++c_ctr;
    }

    for (int c = 0 ; c < LenUnits ; ++c)
    {
        c_vec.push_back(node_coord_units[c]); ++c_ctr;
    }

    for (int v = 0 ; v < NumVerts ; ++v)
    {
        d_vec.push_back(vertex_vector[v].get_x()); ++d_ctr;
        d_vec.push_back(vertex_vector[v].get_y()); ++d_ctr;
        d_vec.push_back(vertex_vector[v].get_z()); ++d_ctr;
    }

    for (int c = 0 ; c < layout.num_cells() ; ++c)
        for (int f = 0 ; f < layout.num_faces(c+1) ; ++f)
        {
            i_vec.push_back(layout(c+1,f+1)); ++i_ctr;
        }

    for (int c = 0 ; c < NumCells ; ++c)
        for (int v = 0 ; v < FOUR ; ++v)
        {
            i_vec.push_back(cells_vertices[c][v]); ++i_ctr;
        }

    for (int s = 0 ; s < NumSides ; ++s)
        for (int v = 0 ; v < THREE ; ++v)
        {
            i_vec.push_back(sides_vertices[s][v]); ++i_ctr;
        }

    if (NumNodeSets > 0)
    {
        for (MAP_String_SetInt::const_iterator flag = node_sets.begin() ;
                flag != node_sets.end() ; ++flag)
        {
            int slen = (*flag).first.length();
            int mems = (*flag).second.size();

            for (int c = 0 ; c < slen ; ++c)
            {
                c_vec.push_back((*flag).first[c]); ++c_ctr;
            }

            i_vec.push_back(slen); ++i_ctr;
            i_vec.push_back(mems); ++i_ctr;
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                    i != (*flag).second.end(); ++i)
            {
                i_vec.push_back(*i); ++i_ctr;
            }
        }
    } // End if (NumNodeSets > 0)

    if (NumSideSets > 0)
    {
        for (MAP_String_SetInt::const_iterator flag = side_sets.begin() ;
                flag != side_sets.end() ; ++flag)
        {
            int slen = (*flag).first.length();
            int mems = (*flag).second.size();

            for (int c = 0 ; c < slen ; ++c)
            {
                c_vec.push_back((*flag).first[c]); ++c_ctr;
            }

            i_vec.push_back(slen); ++i_ctr;
            i_vec.push_back(mems); ++i_ctr;
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                    i != (*flag).second.end(); ++i)
            {
                i_vec.push_back(*i); ++i_ctr;
            }
        }
    } // End if (NumSideSets > 0)

    if (NumCellSets > 0)
    {
        for (MAP_String_SetInt::const_iterator flag = cell_sets.begin() ;
                flag != cell_sets.end() ; ++flag)
        {
            int slen = (*flag).first.length();
            int mems = (*flag).second.size();

            for (int c = 0 ; c < slen ; ++c)
            {
                c_vec.push_back((*flag).first[c]); ++c_ctr;
            }

            i_vec.push_back(slen); ++i_ctr;
            i_vec.push_back(mems); ++i_ctr;
            for (SetInt::const_iterator i = (*flag).second.begin() ;
                    i != (*flag).second.end(); ++i)
            {
                i_vec.push_back(*i); ++i_ctr;
            }
        }
    } // End if (NumCellSets > 0)

    // A terminating null is not necessary, but is harmless,
    // and could improve client safety.
    c_vec.push_back('\0'); ++c_ctr;

    double *d_arr = new double[d_ctr];
    for (int d = 0 ; d < d_ctr ; ++d)
        d_arr[d] = d_vec[d];

    int *i_arr = new int[i_ctr];
    for (int i = 0 ; i < i_ctr ; ++i)
        i_arr[i] = i_vec[i];

    char *c_arr = new char[c_ctr];
    for (int c = 0 ; c < c_ctr ; ++c)
        c_arr[c] = c_vec[c];

    return TET_Mesh::Pack(d_ctr, d_arr, i_ctr, i_arr, c_ctr, c_arr);

}   // end TET_Mesh::pack()

}   // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of mc/TET_Mesh.cc
//---------------------------------------------------------------------------//
