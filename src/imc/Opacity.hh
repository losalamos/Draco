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
    typename MT::CCSF_double sigma_abs;

  // Plankian opacity
    typename MT::CCSF_double planck;

  // fleck factors
    typename MT::CCSF_double fleck;
    
public:
  // opacity constructor
    inline Opacity(const typename MT::CCSF_double &,
		   const typename MT::CCSF_double &,
		   const typename MT::CCSF_double &);

  // member set and accessor functions

    double get_sigma_abs(int cell) const { return sigma_abs(cell); }
    double get_planck(int cell) const { return planck(cell); }
    double get_fleck(int cell) const { return fleck(cell); }
    int num_cells() const { return sigma_abs.get_Mesh().num_cells(); }

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
// constructor

template<class MT>
inline Opacity<MT>::Opacity(const typename MT::CCSF_double &sigma_abs_,
			    const typename MT::CCSF_double &planck_,
			    const typename MT::CCSF_double &fleck_)
    : sigma_abs(sigma_abs_), planck(planck_), fleck(fleck_)
{}

//---------------------------------------------------------------------------//
// Fleck effective scatter cross section

template<class MT>
inline double Opacity<MT>::get_sigeffscat(int cell) const 
{
    return (1.0 -fleck(cell)) * sigma_abs(cell);
}

//---------------------------------------------------------------------------//
// Fleck effective absorption cross section

template<class MT>
inline double Opacity<MT>::get_sigeffabs(int cell) const 
{ 
    return fleck(cell) * sigma_abs(cell);
}

CSPACE

#endif                          // __imc_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.hh
//---------------------------------------------------------------------------//
