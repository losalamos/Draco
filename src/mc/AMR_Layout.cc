//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/AMR_Layout.cc
 * \author Thomas M. Evans
 * \date   Tue Jul 18 16:07:35 2000
 * \brief  AMR_Layout member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "AMR_Layout.hh"
#include <iomanip>

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS 
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.

 * This constructor allows the user to set up the data inside of AMR_Layout
 * based on the number of cells in the layout. The number of faces per cell
 * and the number of cells across a face are as yet undetermined.

 * \param num_cells[=0] number of cells in the problem, defaults to zero

 */
AMR_Layout::AMR_Layout(int num_cells)
    : cell_face(num_cells)
{
    Ensure (cell_face.size() == num_cells);
}

//---------------------------------------------------------------------------//

/*!
 * \brief Constructor when the problem dimension is known up front.

 * This constructor allows the user to enter the number of coarse faces in
 * each cell of the mesh for which this layout is designed.  The internal
 * data is then appropriately constructed based on the size of the mesh and
 * the number of coarse faces per cell.  A drawback when using this
 * constructor is that each cell must have the same number of coarse faces.
 * However, this is not a limitation of the layout. But, if each cell has a
 * different number of coarse faces then they will have to be sized using the
 * set_faces service.

 * Each coarse face is given one fine face by default.

 * \param num_cells number of cells in the layout

 * \param num_faces number of coarse faces per cell

 */
AMR_Layout::AMR_Layout(int num_cells, int num_faces)
    : cell_face(num_cells, vf_int(num_faces, sf_int(1)))
{
    Ensure (cell_face.size() == num_cells);
    Ensure (cell_face[num_cells-1].size() == num_faces);
    Ensure (cell_face[num_cells-1][num_faces-1].size() == 1);
}

//---------------------------------------------------------------------------//
// PUBLIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out the AMR_Layout.
 * \param cell print out layout for given cell index
 */
void AMR_Layout::print(std::ostream &out, int cell) const
{
    using std::endl;
    using std::setw;

    Require (cell > 0 && cell <= cell_face.size());

    int faces = cell_face[cell-1].size();

    // print out the cell info
    out << "====================" << endl;
    out << "Cell         : " << cell << endl;
    if (faces)
	out << "Coarse Faces : " << faces << endl;
    out << "--------------------" << endl;
    for (int i = 0; i < faces; i++)
    {
	out << setw(4) << i+1 << endl;
	for (int j = 0; j < cell_face[cell-1][i].size(); j++)
	    out << setw(8) << j+1 << setw(5) 
		<< cell_face[cell-1][i][j] << endl;
    }
    out << "====================" << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test for equality.
 */
bool AMR_Layout::operator==(const AMR_Layout &rhs) const
{
    // if the cell data is equal then return true
    if (cell_face == rhs.cell_face)
	return true;

    return false;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out all cells using operator<<.
 */
std::ostream& operator<<(std::ostream &out, const AMR_Layout &object)
{
    int num_cells = object.num_cells();
    for (int i = 1; i <= num_cells; i++)
	object.print(out, i);
    return out;
}

} // end namespace rtt_mc


//---------------------------------------------------------------------------//
//                              end of AMR_Layout.cc
//---------------------------------------------------------------------------//
