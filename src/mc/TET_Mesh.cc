//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/TET_Mesh.cc
 * \author H. Grady Hughes
 * \date   Tue Jan 18 10:15:43 MST 2000
 * \brief  Tetrahedral mesh class implementation file.
 */
//---------------------------------------------------------------------------//

#include "TET_Mesh.hh"
#include "Constants.hh"
#include <iomanip>

namespace rtt_mc
{

//___________________________________________________________________________//
/*!
 * \brief       TET_Mesh constructor.
 * \param coord_          Smart pointer to base class of derived XYZ system.
 * \param layout_         Layout of mesh.
 * \param vertex_vector_  ThreeVector vertices in a list of all vertices.
 * \param cells_vertices_ Internal identifiers of the four vertices of cells.
 * \param submesh_        Submesh indicator flag.
 */
TET_Mesh::TET_Mesh(SP<Coord_sys> coord_, Layout & layout_,
    SF_THREEVECTOR & vertex_vector_, VF_INT & cells_vertices_, bool submesh_)
    : coord(coord_), layout(layout_), vertex_vector(vertex_vector_),
      cells_vertices(cells_vertices_), submesh(submesh_)
{
    // For a TET_Mesh, there must be an XYZ coordinate system.
    Check (coord);
    Check (THREE == coord->get_dim());
    Check (string("xyz") == coord->get_Coord());

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
 * \brief          Determine whether a given position is inside a given cell.
 * \param position XYZ-position as STL vector.
 * \param cell     External number of cell.
 * \return         True if position is strictly within the cell, else false.
 *
 * The external number "cell" is checked for validity, then converted to
 * cell_ = cell - 1 for internal use.
 */
bool TET_Mesh::in_cell(const SF_DOUBLE &position, int cell) const
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

}   // end TET_Mesh::in_cell(const SF_DOUBLE &, int)

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
double TET_Mesh::get_db(const SF_DOUBLE &position, const SF_DOUBLE &omega,
                        int cell, int &face) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);
    ThreeVector UVW(omega);
    SF_DOUBLE dist(FOUR);

    for (int f_ = 0 ; f_ < FOUR ; f_++)
        {
            ThreeVector N = get_outward_cross(cell_, f_);
            double denom = N.dot(UVW);
            int v = (f_ + 1) % FOUR;  // any vertex on the face.

            if (denom > 0.0)
                dist.push_back(N.dot(
                    vertex_vector[cells_vertices[cell_][v]] - XYZ)/denom);
            else
                dist.push_back(global::huge);
        }

    SF_DOUBLE::iterator itor = min_element(dist.begin(),dist.end());

    double dist_min = *itor;
    Ensure ( dist_min > 0.0 && dist_min < global::huge );

    int face_ = itor - dist.begin();
    face = face_ + 1;
    Ensure ( face >= 1 && face <= FOUR );

    return *itor;

}   // end TET_Mesh::get_db(const SF_DOUBLE&,const SF_DOUBLE&,int,int&)

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
const TET_Mesh::SF_DOUBLE TET_Mesh::get_normal(int cell, int face) const
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
const TET_Mesh::SF_DOUBLE TET_Mesh::get_normal_in(int cell, int face) const
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
 * dim# = 0, 1, 2 <====> for axis X, Y, Z.
 *
 * vertex# = 0, 1, 2, 3.
 */
const TET_Mesh::VF_DOUBLE TET_Mesh::get_vertices(int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;
    // Could check that (cells_vertices[cell_].size() == FOUR);

    TET_Mesh::VF_DOUBLE ret_vert(THREE);

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
 * dim# = 0, 1, 2 <====> for axis X, Y, Z.
 *
 * vertex# = 0, 1, 2, but representing three out of four of (0, 1, 2, 3),
 * since each face has 3 vertices.  For no particular reason, these are
 * chosen so that right-handed circulation in the order specified would
 * result in an outward direction.
 */
const TET_Mesh::VF_DOUBLE TET_Mesh::get_vertices(int cell, int face) const
{
    Valid(cell, face);
    int cell_ = cell - 1;
    int face_ = face - 1;

    VF_DOUBLE ret_vert(THREE);

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
 */
const TET_Mesh::SF_DOUBLE TET_Mesh::sample_pos(int cell, Sprng &random) const
{
    Valid(cell);
    int cell_ = cell - 1;

    double fraction = std::pow(random.ran(),1.0/3.0);
    Check ( fraction > 0.0 && fraction < 1.0 );

    int v0 = cells_vertices[cell_][0];
    int v1 = cells_vertices[cell_][1];
    int v2 = cells_vertices[cell_][2];
    int v3 = cells_vertices[cell_][3];

    ThreeVector A = lin_comb(vertex_vector[v3],vertex_vector[v0],fraction);
    ThreeVector B = lin_comb(vertex_vector[v3],vertex_vector[v1],fraction);
    ThreeVector C = lin_comb(vertex_vector[v3],vertex_vector[v2],fraction);

    return sample_in_triangle(A, B, C, random).convert();

}   // end TET_Mesh::sample_pos(int,Sprng &)

//___________________________________________________________________________//
/*!
 * \brief           Sample position in a given cell with a tilt.
 * \param cell      External number of cell.
 * \param random    The random-number generator, used as random.ran()
 * \param slope     Specifies a tilt or slope-like function.
 * \param center_pt Reference location for tilt or slope-like function.
 * \return          Scalar_field[dim#] == sampled coordinate along dim#-axis.
 */
const TET_Mesh::SF_DOUBLE TET_Mesh::sample_pos(int cell, Sprng &random,
    SF_DOUBLE slope, double center_pt) const
{
    // QUESTION: Should "slope" be "const SF_DOUBLE &slope" ?
    // IMPLEMENTATION LATER.  For the moment, ignore the tilt.

    return TET_Mesh::sample_pos(cell, random);  // Uniform sampling.

}   // end TET_Mesh::sample_pos(int,Sprng &,SF_DOUBLE,double)

//___________________________________________________________________________//
/*!
 * \brief        Sample position uniformly on a given face of a given cell.
 * \param cell   External number of cell.
 * \param face   External identifier for the face (1, 2, 3, or 4).
 * \param random The random-number generator, used as random.ran()
 * \return       Scalar_field[dim#] == sampled coordinate along dim#-axis.
 */
const TET_Mesh::SF_DOUBLE TET_Mesh::sample_pos_on_face(int cell, int face,
    Sprng &random) const
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

    return sample_in_triangle(vertex_vector[cells_vertices[cell_][va]],
                              vertex_vector[cells_vertices[cell_][vb]],
                              vertex_vector[cells_vertices[cell_][vc]],
                              random).convert();

}   // end TET_Mesh::sample_pos_on_face(int,int,Sprng &)

//___________________________________________________________________________//
/*!
 * \brief Overloaded operator== for design-by-contract.
 */
bool TET_Mesh::operator==(const TET_Mesh &rhs) const
{
    // Verify that we have the same coordinate systems.
//  (After the coordinate system is installed.)
//  if (coord != rhs.coord)
//  return false;

// QUESTION: Should we check submesh ?

    // Verify that the Layouts are equal.
    if (layout != rhs.layout)
        return false;

    // Verify the coordinates of the vertices.
    if (vertex_vector != rhs.vertex_vector)
        return false;

    // Verify the identities of the vertices of each cell.
    if (cells_vertices != rhs.cells_vertices)
        return false;

    // if we haven't returned, then the two meshes must be equal
    return true;
}   // end TET_Mesh::operator==(const TET_Mesh &)

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
int TET_Mesh::get_cell(const SF_DOUBLE &position) const
{
    Require (!submesh);

    int ret_cell = 0;

    for (int cell = 1 ; cell <= num_cells() ; cell++)
        if ( in_cell(position, cell) )
            ret_cell = cell;

    Ensure (ret_cell > 0);
    return ret_cell;

}   // end TET_Mesh::get_cell(const SF_DOUBLE &)

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
 */
double TET_Mesh::get_min_db(const SF_DOUBLE &position, int cell) const
{
    Valid(cell);
    int cell_ = cell - 1;

    ThreeVector XYZ(position);
    SF_DOUBLE dist(FOUR);

    for (int face = 1 ; face <= FOUR ; face++)
        {
            ThreeVector N(get_normal(cell, face));
            int v_ = face % FOUR;  // any vertex on the face.

            dist.push_back(
                N.dot(vertex_vector[cells_vertices[cell_][v_]] - XYZ));
        }

    SF_DOUBLE::iterator itor = min_element(dist.begin(),dist.end());

    double dist_min = *itor;
    Ensure ( dist_min > 0.0 && dist_min < global::huge );
    return dist_min;

}   // end TET_Mesh::get_min_db(const SF_DOUBLE &, int)

}   // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of mc/TET_Mesh.cc
//---------------------------------------------------------------------------//
