//----------------------------------*-C++-*----------------------------------//
// Layout.cc
// Thomas M. Evans
// Fri Jan 30 15:53:52 1998
//---------------------------------------------------------------------------//
// @> Layout class implementation file
//---------------------------------------------------------------------------//

#include "Layout.hh"
#include <iomanip>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// print member function used for printing the cell-face-cell info.
// for one particular cell
void Layout::Print(int cell_index) const
{
    using std::cout;
    using std::endl;
    using std::setw;
      
    int num_faces = Num_faces(cell_index);

  // print the face information for cell indexed by cell_index
    cout << "=============" << endl;
    cout << "Cell   : "     << cell_index << endl;
    cout << "# faces: "     << num_faces << endl;
    cout << "-------------" << endl;
    cout << setw(4) << "Face" << setw(9) << " Neighbor" << endl;
    cout << "-------------" << endl;
    for (int i = 0; i <= num_faces-1; i++)
        cout << setw(4) << i+1 << setw(9)
             << face_cell[cell_index-1][i] << endl;
    cout << "=============" << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
ostream & operator<<(ostream &output, const Layout &object)
{
    int num_cells = object.Num_cells();
    for (int i = 1; i <= num_cells; i++)
        object.Print(i);
    return output;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Layout.cc
//---------------------------------------------------------------------------//
