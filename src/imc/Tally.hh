//----------------------------------*-C++-*----------------------------------//
// Tally.hh
// Todd J. Urbatsch
// Mon Apr  6 14:38:03 1998
//---------------------------------------------------------------------------//
// @> Tally class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Tally_hh__
#define __imc_Tally_hh__

//===========================================================================//
// class Tally - 
//
// Purpose : Accumulate Particle tallies.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "ds++/SP.hh"
#include <iostream>

IMCSPACE

using std::ostream;

template<class MT>
class Tally 
{

private:
    typename MT::CCSF_double energy_dep;
    double energy_dep_tot;

    typename MT::CCSF_double eweighted_pathlen;
    typename MT::CCSF_double census_energy;
    double new_ecen_tot;

public:
  // Tally constructor
    explicit Tally(SP<MT> mesh) 
	: energy_dep(mesh), energy_dep_tot(0.0), eweighted_pathlen(mesh), 
	  census_energy(mesh), new_ecen_tot(0.0){}

  // accumulate energy deposited
    void deposit_energy(const int cell, const double energy);

  // access energy deposited for a cell, total.
    inline double get_energy_dep(const int cell);
    inline double get_energy_dep_tot();

  // accumulate energy-weighted path length
    void accumulate_ewpl(const int cell, const double ewpl);

  // accumulate new census energy
    void accumulate_ecen(const int cell, const double new_ecen);

  // access accumulated energy-weighted path length per cell
    inline double get_accum_ewpl(const int cell);

  // access accumulated new census energy per cell, total
    inline double get_new_ecen(const int cell);
    inline double get_new_ecen_tot();


  // accessors
    int num_cells() const { return energy_dep.get_Mesh().num_cells(); }

  // diagnostics for tally
    void print(ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &out, const Tally<MT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// inline member functions for Tally
//---------------------------------------------------------------------------//

template<class MT>
inline double Tally<MT>::get_energy_dep(const int cell)
{
    return energy_dep(cell);
}

//---------------------------------------------------------------------------//

template<class MT>
inline double Tally<MT>::get_energy_dep_tot()
{
    return energy_dep_tot;
}

//---------------------------------------------------------------------------//

template<class MT>
inline double Tally<MT>::get_accum_ewpl(const int cell)
{
    return eweighted_pathlen(cell);
}

//---------------------------------------------------------------------------//

template<class MT>
inline double Tally<MT>::get_new_ecen(const int cell)
{
    return census_energy(cell);
}

//---------------------------------------------------------------------------//

template<class MT>
inline double Tally<MT>::get_new_ecen_tot()
{
    return new_ecen_tot;
}


CSPACE

#endif                          // __imc_Tally_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Tally.hh
//---------------------------------------------------------------------------//

