//----------------------------------*-C++-*----------------------------------//
// Opacity.hh
// Thomas M. Evans
// Fri Feb  6 13:52:29 1998
//---------------------------------------------------------------------------//
// @> Opacity class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Opacity_hh__
#define __imc_Opacity_hh__

//===========================================================================//
// class Opacity - 
//
// Date created : 2-6-98
// Purpose      : Simple opacity package for IMCTEST
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
class Opacity
{
private:
  // sigma = kappa * rho
    typename MT::CCSF_double sigma;

  // Plankian opacity
    typename MT::CCSF_double planck;

  // fleck factors
    typename MT::CCSF_double fleck;
    
public:
  // opacity constructor
    explicit Opacity(const typename MT::CCSF_double &sigma_) 
	: sigma(sigma_), planck(sigma_), fleck(sigma_) {}

  // member set and accessor functions

    double get_sigma(int cell) const { return sigma(cell); }
    double get_planck(int cell) const { return planck(cell); }
    double get_fleck(int cell) const { return fleck(cell); }
    int num_cells() const { return sigma.get_Mesh().num_cells(); }

  // operations
    double fplanck(int cell) const { return fleck(cell) * planck(cell); }
    inline double get_sigeffscat(int cell) const;
    inline double get_sigeffabs(int cell) const;

  // diagnostic member functions
    void print(ostream &, int) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &, const Opacity<MT> &);

//---------------------------------------------------------------------------//
// inline member functions for Opacity
//---------------------------------------------------------------------------//

template<class MT>
inline double Opacity<MT>::get_sigeffscat(int cell) const 
{
    return fleck(cell) * sigma(cell);
}

//---------------------------------------------------------------------------//

template<class MT>
inline double Opacity<MT>::get_sigeffabs(int cell) const 
{ 
    return (1.0 - fleck(cell)) * sigma(cell);
}

CSPACE

#endif                          // __imc_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.hh
//---------------------------------------------------------------------------//
