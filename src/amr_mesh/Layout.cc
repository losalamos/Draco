//----------------------------------*-C++-*----------------------------------//
// Layout.cc
// Todd Adams
// Thu 23 Sept 15:53:52 1998
/*! 
 * \file   amr_mesh/Layout.cc
 * \author B.T. Adams
 * \date   Wed Sep 15 10:33:26 1999
 * \brief  Implementation file for Layout class library.
 */
//---------------------------------------------------------------------------//
// @> Layout class implementation file - this file should really be moved 
// back to mc. It was modified from the original for AMR meshes to allow a
// cell to have more than one neighbor, but it defaults to a single adjacent
// cell and thus should work just hunky dory with the OS_Mesh type too (with
// a small change to OS_Builder.cc also included to accomodate a ragged right
// array for the adjacent cell vector).
//---------------------------------------------------------------------------//

#include "Layout.hh"
#include <iomanip>

namespace rtt_amr 
{

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// print member function used for printing the cell-face-cell info.
// for one particular cell
void Layout::print(ostream & output, int cell_index) const
{
    using std::endl;
    using std::setw;
      
    int faces = get_num_cell_faces(cell_index);

    // print the adjacent cell information for face indexed by face_index of 
    // cell indexed by cell_index
    output << "=============" << endl;
    output << "Cell   : "     << cell_index << endl;
    output << "# faces: "     << faces << endl;
    output << "-------------" << endl;
    output << setw(4) << "Face" << setw(9) << " Neighbor" << endl;
    output << "-------------" << endl;
    for (int face_index = 0; face_index < faces; face_index++)
    {
        output << setw(4) << face_index;
        int adjCells = get_num_adj_cells(cell_index,face_index);
	for (int adj = 0; adj < adjCells; adj++)
	    output << setw(9) << face_cell[cell_index][face_index][adj] 
		   << endl;
    }
    output << "=============" << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
ostream & operator<<(ostream & output, const Layout & object)
{
    int num_cells = object.get_num_cells();
    for (int cell_index = 0; cell_index < num_cells; cell_index++)
        object.print(output, cell_index);
    return output;
}

} // end namespace rtt_amr

//---------------------------------------------------------------------------//
//                              end of Layout.cc
//---------------------------------------------------------------------------//
