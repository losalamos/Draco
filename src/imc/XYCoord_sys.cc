//----------------------------------*-C++-*----------------------------------//
// XYCoord_sys.cc
// Thomas M. Evans
// Fri Jan 30 16:52:13 1998
//---------------------------------------------------------------------------//
// @> XYCoord_sys derived class implementation file
//---------------------------------------------------------------------------//

#include "imctest/XYCoord_sys.hh"
#include "imctest/Global.hh"
#include "imctest/Random.hh"
#include <iostream>
#include <cmath>

IMCSPACE

//---------------------------------------------------------------------------//
// virtual member functions
//---------------------------------------------------------------------------//
// set Omega directions in a 2D XY geometry
void XYCoord_sys::Set_omega(vector<double> &omega_, Random &random) const
{
  // sample phi for 2D XY coordinate system
    double costheta;
    costheta  = 1.0 - 2.0 * random.ran();

  // calculate 2D direction cosines
    omega_[0] = costheta;
    omega_[1] = sqrt(1.0-costheta*costheta);
}

// calculate Omega directions in 2D XY geometry after a scatter
void XYCoord_sys::Calc_omega(double costheta, double phi, vector<double> 
			     &omega_) const
{
    using Global::pi;
    using std::cos;
    using std::sin;

  // calculate new direction cosines
    double theta     = acos(costheta);
    double new_theta = theta + acos(omega_[0]);
    omega_[0]        = cos(new_theta);
    omega_[1]        = sin(new_theta);
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of XYCoord_sys.cc
//---------------------------------------------------------------------------//
