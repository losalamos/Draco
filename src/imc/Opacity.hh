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
  // cell-centered array opacities
    typename MT::CCSF sigma;
public:
  // opacity constructor
    explicit Opacity(const typename MT::CCSF &sigma_) : sigma(sigma_) {}

  // member set and accessor functions
    double& Sigma(int cell) { return sigma(cell); }
    double Sigma(int cell) const { return sigma(cell); }
    int Num_cells() const { return sigma.Mesh().Num_cells(); }

  // diagnostic member functions
    void Print(ostream &, int) const;
};

// overloaded operators
template<class MT>
ostream& operator<<(ostream &, const Opacity<MT> &);

CSPACE

#endif                          // __imctest_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Opacity.hh
//---------------------------------------------------------------------------//
