//----------------------------------*-C++-*----------------------------------//
// Opacity.hh
// Thomas M. Evans
// Fri Feb  6 13:52:29 1998
//---------------------------------------------------------------------------//
// @> Opacity class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Opacity_hh__
#define __imctest_Opacity_hh__

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

#include "imctest/Names.hh"
#include "ds++/SP.hh"
#include <iostream>

IMCSPACE

using std::ostream;

template<class MT>
class Opacity
{
private:
  // sigma = kappa * rho
    typename MT::CCSF<double> sigma;

  // Plankian opacity
    typename MT::CCSF<double> planck;

  // fleck factors
    typename MT::CCSF<double> fleck;
    
public:
  // opacity constructor
    explicit Opacity(const typename MT::CCSF<double> &sigma_) 
	: sigma(sigma_), planck(sigma_), fleck(sigma_) {}

  // member set and accessor functions

    double& get_sigma(int cell) { return sigma(cell); }
    double get_sigma(int cell) const { return sigma(cell); }
    int num_cells() const { return sigma.get_Mesh().num_cells(); }

  // operations
    double fplanck(int cell) { return fleck(cell) * planck(cell); }

  // diagnostic member functions
    void print(ostream &, int) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
ostream& operator<<(ostream &, const Opacity<MT> &);

CSPACE

#endif                          // __imctest_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Opacity.hh
//---------------------------------------------------------------------------//
