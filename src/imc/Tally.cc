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
    energy_dep_tot   += energy;
}

template<class MT>
void Tally<MT>::accumulate_ewpl(const int cell, const double ewpl)
{
    eweighted_pathlen(cell) += ewpl;
}


template<class MT>
void Tally<MT>::accumulate_cen_info(const int cell, const double new_ecen)
{
    census_energy(cell) += new_ecen;
    new_ecen_tot        += new_ecen;
    new_ncen(cell)      += 1;
    new_ncen_tot        += 1;
}

//---------------------------------------------------------------------------//
// print functions
//---------------------------------------------------------------------------//

template<class MT>
void Tally<MT>::print(ostream &out) const
{
    out << setw(8) << setiosflags(ios::right) << "Cell"
	<< setw(15) << "energy-dep" << setw(15) << "ewpl" 
	<< setw(15) << "new_ecen" 
	<< setw(15) << "new_ncen" << endl;
    out << "------------------------------------" <<
           "------------------------------------" << endl;

    out.precision(4);
    for (int i = 1; i <= energy_dep.get_Mesh().num_cells(); i++)
	out << setw(8) << i 
	    << setw(15) << setiosflags(ios::scientific) << energy_dep(i) 
	    << setw(15) << setiosflags(ios::scientific) << eweighted_pathlen(i) 
	    << setw(15) << setiosflags(ios::scientific) << census_energy(i) 
	    << setw(15) << setiosflags(ios::scientific) << new_ncen(i) 
	    << endl;

    out << endl;
    out << "Total new census energy: " 
	<< setw(15) << setiosflags(ios::scientific) << new_ecen_tot
	<< ", new_ncen_tot: " << new_ncen_tot << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Tally.cc
//---------------------------------------------------------------------------//
