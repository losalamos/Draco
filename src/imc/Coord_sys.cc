//----------------------------------*-C++-*----------------------------------//
// Coord_sys.cc
// Thomas M. Evans
// Fri Jan 30 16:45:37 1998
//---------------------------------------------------------------------------//
// @> Coord_sys base class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Coord_sys.hh"
#include "imctest/Global.hh"
#include "imctest/Random.hh"
#include <iostream>
#include <cmath>

IMCSPACE

//---------------------------------------------------------------------------//
// virtual member functions
//---------------------------------------------------------------------------//
// set Omega directions for 3D transport
void Coord_sys::Set_omega(vector<double> &omega_, Random &random) const
{
    using Global::pi;
    using std::cos;
    using std::sin;
    using std::sqrt;

  // sample costheta and phi for 3D transport 
    double costheta, sintheta, phi;
    costheta = 1 - 2 * random.ran();
    sintheta = sqrt(1 - costheta * costheta);
    phi      = 2 * pi * random.ran();

  // calculate 3D direction cosines
    omega_[0] = sintheta * cos(phi);
    omega_[1] = sintheta * sin(phi);
    omega_[2] = costheta;
}

// calculate Omega directions for 3D transport
void Coord_sys::Calc_omega(double costheta, double phi, vector<double> 
			      &omega_) const
{
    using Global::pi;
    using std::cos;
    using std::sin;
    using std::sqrt;

  // calculate new direction cosines
    double sintheta = sqrt(1 - costheta * costheta);
    vector<double> old_dir;
    old_dir = omega_;
    double factor = sqrt(1 - old_dir[2] * old_dir[2]);
    if (factor == 0)
    {
        omega_[0] = sintheta * cos(phi);
        omega_[1] = sintheta * sin(phi);
        omega_[2] = old_dir[2] * cos(phi);
    }
    else
    {
        omega_[0] = old_dir[0] * costheta + old_dir[2] * old_dir[0] *
            sintheta *cos(phi) / factor - old_dir[1] * sintheta *
            sin(phi) / factor;
        omega_[1] = old_dir[1] * costheta + old_dir[2] * old_dir[1] *
            sintheta *cos(phi) / factor + old_dir[0] * sintheta *
            sin(phi) / factor;
        omega_[2] = old_dir[2] * costheta - factor * sintheta * cos(phi);
    }
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Coord_sys.cc
//---------------------------------------------------------------------------//
