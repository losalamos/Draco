//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Coord_sys.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 30 16:45:37 1998
 * \brief  Coord_sys base class implementation file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Coord_sys.hh"
#include "Constants.hh"
#include "Math.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_mc 
{

using global::soft_equiv;
using global::dot;
using global::pi;
using std::cos;
using std::sin;
using std::sqrt;
using std::vector;
using std::string;

//---------------------------------------------------------------------------//
// VIRTUAL MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
// set Omega directions for 3D transport

Coord_sys::sf_double Coord_sys::sample_dir(std_string dist, 
					   rng_Sprng &random) const
{
    // make return vector
    vector<double> omega_(3);

    // get direction cosines for different distributions
    if (dist == "isotropic")
    {
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

    Check (soft_equiv(dot(omega_, omega_), 1.0, 1.e-6));

    // return vector
    return omega_;
}

//---------------------------------------------------------------------------//
// calculate Omega directions for 3D transport

void Coord_sys::calc_omega(double costheta, double phi, 
			   sf_double &omega_) const
{
    // calculate new direction cosines
    double sintheta = sqrt(1 - costheta * costheta);
    vector<double> old_dir;
    old_dir = omega_;
    double factor = sqrt(1 - old_dir[2] * old_dir[2]);
    if (factor < 1.e-6)
    {
        omega_[0] = sintheta * cos(phi);
        omega_[1] = sintheta * sin(phi);
        omega_[2] = old_dir[2] / std::fabs(old_dir[2]) * costheta;
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

    Check (soft_equiv(dot(omega_, omega_), 1.0, 1.e-6));
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Coord_sys.cc
//---------------------------------------------------------------------------//
