//----------------------------------*-C++-*----------------------------------//
// Tally.cc
// Todd J. Urbatsch
// Mon Apr  6 14:38:03 1998
//---------------------------------------------------------------------------//
// @> Tally class implementation file.
//---------------------------------------------------------------------------//

#include "imc/Tally.hh"
#include <iomanip>

IMCSPACE

using std::setw;
using std::ios;
using std::setiosflags;

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//

template<class MT>
void Tally<MT>::deposit_energy(const int cell, const double energy)
{
    energy_dep(cell) += energy;
    energy_dep_tot += energy;
}

//---------------------------------------------------------------------------//
// print functions
//---------------------------------------------------------------------------//

template<class MT>
void Tally<MT>::print(ostream &out) const
{
    out << setw(8) << setiosflags(ios::right) << "Cell"
	<< setw(15) << "Tally" << endl;
    out << "-----------------------" << endl;

    out.precision(4);
    for (int i = 1; i <= energy_dep.get_Mesh().num_cells(); i++)
	out << setw(8) << i << setw(15) << setiosflags(ios::scientific)
	    << energy_dep(i) << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Tally.cc
//---------------------------------------------------------------------------//
