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

#include "Names.hh"
#include "SP.hh"

IMCSPACE

template<class MT>
class Opacity
{
private:
    typename MT::CCSF sigma;
    typename MT::CCSF sigma_a;
    typename MT::CCSF sigma_s;
public:
    Opacity(SP<MT> mesh) : sigma(mesh), sigma_a(mesh), sigma_s(mesh)
    {}
    void Quickset(double, double, double);
  // member set and accessor functions
    double& Sigma(int cell) { return sigma(cell); }
    double Sigma(int cell) const { return sigma(cell); }
    double& Sigma_a(int cell) { return sigma_a(cell); }
    double Sigma_a(int cell) const { return sigma_a(cell); }
    double& Sigma_s(int cell) { return sigma_s(cell); }
    double Sigma_s(int cell) const { return sigma_s(cell); }
};

CSPACE

#endif                          // __imctest_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Opacity.hh
//---------------------------------------------------------------------------//
