//----------------------------------*-C++-*----------------------------------//
// Tally.hh
// Todd J. Urbatsch
// Mon Apr  6 14:38:03 1998
//---------------------------------------------------------------------------//
// @> Tally class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Tally_hh__
#define __imctest_Tally_hh__

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

#include "imctest/Names.hh"
#include "ds++/SP.hh"



IMCSPACE

template<class MT>
class Tally 
{

private:
    typename MT::CCSF_double energy_dep;
    double energy_dep_tot;

public:
  // Tally constructor
    explicit Tally(SP<MT> mesh) 
	: energy_dep(mesh), energy_dep_tot(0.0){}

  // accumulate energy deposited
    void deposit_energy(const int cell, const double energy);

  // access energy deposited for a cell, total.
    inline double get_energy_dep(const int cell);
    inline double get_energy_dep_tot();
};

//---------------------------------------------------------------------------//
// inline member functions for Tally
//---------------------------------------------------------------------------//

 
template<class MT>
inline double Tally<MT>::get_energy_dep(const int cell)
{
    return energy_dep(cell);
}

template<class MT>
inline double Tally<MT>::get_energy_dep_tot()
{
    return energy_dep_tot;
}
  


CSPACE

#endif                          // __imctest_Tally_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Tally.hh
//---------------------------------------------------------------------------//

