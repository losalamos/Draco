//----------------------------------*-C++-*----------------------------------//
// Mat_State.cc
// Thomas M. Evans
// Mon Mar  9 16:06:28 1998
//---------------------------------------------------------------------------//
// @> Mat_State class implementation file
//---------------------------------------------------------------------------//

#include "imc/Mat_State.hh"
#include <iomanip>

IMCSPACE

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//

template<class MT>
void Mat_State<MT>::print(ostream &output, int cell) const
{
  // print out material state of cell
    using std::cout;
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output.precision(4);
       
    output << setw(8) << cell << setw(15) << setiosflags(ios::scientific)
	   << density(cell) << setw(15) << temperature(cell) << endl;
} 
//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &output, const Mat_State<MT> &object)
{
  // print out opacities for all cells
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;
    
    output << setw(8)  << setiosflags(ios::right) << "Cell" 
	   << setw(15) << "Density" << setw(15) << "Temp" << endl;
    output << "--------------------------------------" << endl;

    for (int i = 1; i <= object.num_cells(); i++)
	object.print(output, i);
    return output;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Mat_State.cc
//---------------------------------------------------------------------------//
