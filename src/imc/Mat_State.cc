//----------------------------------*-C++-*----------------------------------//
// Mat_State.cc
// Thomas M. Evans
// Mon Mar  9 16:06:28 1998
//---------------------------------------------------------------------------//
// @> Mat_State class implementation file
//---------------------------------------------------------------------------//

#include "Mat_State.hh"
#include <iomanip>

IMCSPACE

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
template<class MT>
void Mat_State<MT>::Print(int cell) const
{
  // print out material state of cell
    using std::cout;
    using std::endl;
    using std::setw;

    cout << setw(8) << cell << setw(10) << density(cell) << setw(10) 
	 << temperature(cell) << endl;
} 

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
template<class MT>
ostream& operator<<(ostream &output, const Mat_State<MT> &object)
{
  // print out opacities for all cells
    using std::cout;
    using std::endl;
    using std::setw;

    cout << "  Cell  " << " Density  " << "   Temp   " << endl;
    cout << "----------------------------------------" << endl;

    for (int i = 1; i <= Num_cells(); i++)
        object.Print(i);
    return output;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Mat_State.cc
//---------------------------------------------------------------------------//
