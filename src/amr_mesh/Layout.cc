//----------------------------------*-C++-*----------------------------------//
// Layout.cc
// Thomas M. Evans
// Fri Jan 30 15:53:52 1998
//---------------------------------------------------------------------------//
// @> Layout class implementation file
//---------------------------------------------------------------------------//

#include "Layout.hh"
#include <iomanip>

namespace rtt_mc 
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
void Layout::print(ostream &output, int cell_index) const
{
    using std::endl;
    using std::setw;
      
    int faces = num_faces(cell_index);

    // print the adjacent cell information for face indexed by face_index of 
    // cell indexed by cell_index
    output << "=============" << endl;
    output << "Cell   : "     << cell_index << endl;
    output << "# faces: "     << faces << endl;
    output << "-------------" << endl;
    output << setw(4) << "Face" << setw(9) << " Neighbor" << endl;
    output << "-------------" << endl;
    for (int face_index = 1; face_index <= faces; face_index++)
    {
        output << setw(4) << face_index;
        int adjCells = num_adj(cell_index,face_index);
	for (int adj = 1; adj <= adjCells; adj++)
	    output << setw(9) << face_cell[cell_index-1][face_index-1][adj-1] 
		   << endl;
    }
    output << "=============" << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
ostream& operator<<(ostream &output, const Layout &object)
{
    int num_cells = object.num_cells();
    for (int cell_index = 1; cell_index <= num_cells; cell_index++)
        object.print(output, cell_index);
    return output;
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Layout.cc
//---------------------------------------------------------------------------//
