//----------------------------------*-C++-*----------------------------------//
// Opacity.t.hh
// Thomas M. Evans
// Fri Feb  6 13:52:29 1998
//---------------------------------------------------------------------------//
// @> Opacity class implementation file
//---------------------------------------------------------------------------//

#include "Opacity.hh"
#include <iomanip>

namespace rtt_imc 
{

using std::ostream;

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// print function for diagnostics

template<class MT>
void Opacity<MT>::print(ostream &output, int cell) const
{
  // print out opacities in cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output.precision(4);

    output << setw(8) << cell << setw(15) << setiosflags(ios::scientific)
	   << sigma_abs(cell) << setw(15) << sigma_thomson(cell) << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &output, const Opacity<MT> &object)
{
  // print out opacities for all cells
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output << setw(8) << setiosflags(ios::right) << "Cell" 
	   << setw(15) << "Abs Opacity" 
	   << setw(15) << "Thomson Opac" << endl;
    output << "--------------------------------------" << endl;

    for (int i = 1; i <= object.num_cells(); i++)
        object.print(output, i);
    return output;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Opacity.t.hh
//---------------------------------------------------------------------------//
