//----------------------------------*-C++-*----------------------------------//
// Opacity.cc
// Thomas M. Evans
// Fri Feb  6 13:52:29 1998
//---------------------------------------------------------------------------//
// @> Opacity class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Opacity.hh"
#include <iomanip>

IMCSPACE

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// print function for diagnostics
template<class MT>
void Opacity<MT>::Print(int cell) const
{
  // print out opacities in cell
    using std::cout;
    using std::endl;
    using std::setw;

    cout << setw(8) << cell << setw(10) << sigma(cell) << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
template<class MT>
ostream& operator<<(ostream &output, const Opacity<MT> &object)
{
  // print out opacities for all cells
    using std::cout;
    using std::endl;
    using std::setw;

    cout << "  Cell  " << " Opacity  " << endl;
    cout << "------------------------" << endl;

    for (int i = 1; i <= Num_cells(); i++)
        object.Print(i);
    return output;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity.cc
//---------------------------------------------------------------------------//
