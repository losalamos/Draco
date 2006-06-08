//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/Ylm.cc
 * \author Kent Budge
 * \date   Tue Sep 21 09:20:10 2004
 * \brief  Implementation of Ylm
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include "ds++/Assert.hh"
#include "units/PhysicalConstants.hh"
#include "Ylm.hh"

namespace rtt_sf
{
using rtt_units::PI;

//---------------------------------------------------------------------------//
/*!
 * We use the real representation of the spherical harmonics
 *
 * The spherical harmonics \f$ Y_l^m (\theta,\phi) \f$ are the angular portion
 * of the solution to <a
 * href="http://mathworld.wolfram.com/LaplacesEquation.html">Laplace's
 * equation</a> in <a
 * href="http://mathworld.wolfram.com/SphericalCoordinates.html">spherical
 * coordinates</a> where azimuthal symmetry is not present. Some care must be
 * taken in identifying the notational convention being used. In this entry,
 * \f$ \theta \f$ is taken as the polar (colatitudinal) coordinate with \f$
 * \theta\in [0,\pi] \f$, and \f$ \phi \f$ as the azimuthal (longitudinal)
 * coordinate with \f$ \phi\in [0,2\pi ] \f$.
 * 
 * \f[
 * Y_{lm}(\theta,\phi) = (-1)^{\frac{\left| m
 *
 * \right|+m}{2}}i^{m+1} \sqrt{\frac{2l+1}{4\pi}
 * \frac{(l-\left|m\right|)!}{(l+\left|m\right|)!}} P_l^{\left|m\right|} (\cos
 * \theta)\frac{e^{im\phi}-(-1)^me^{-im\phi}}{\sqrt{2}}
 * \qquad for \qquad m<0,
 * \f]
 * \f[
 * Y_{lm}(\theta,\phi) = \sqrt{\frac{2l+1}{4\pi}}P_l^0(\cos \theta) \qquad for
 * \qquad m=0,
 * \f]
 * and
 * \f[
 * Y_{lm}(\theta,\phi) = (-1)^{\frac{\left| m
 * \right|+m}{2}}i^m \sqrt{\frac{2l+1}{4\pi}
 * \frac{(l-m)!}{(l+m)!}} P_l^m (\cos
 * \theta)\frac{e^{im\phi}+(-1)^me^{-im\phi}}{\sqrt{2}} \qquad for \qquad m>0.
 * \f]
 *
 * \pre \f$\left|m\right| \le l\f$
 *
 * \sa http://mathworld.wolfram.com/SphericalHarmonic.html
 */
double Ylm(unsigned const l, int const m, double const theta, double phi)
{
    Require(abs(m)<=l);

    unsigned const absm = abs(m);
    double const cost = cos(theta);
    double Result = gsl_sf_legendre_sphPlm( l, absm, cost ) ;

    // Not sure how important this sign convention really is, but ...
    if (m<0)
    {
        if (absm%2)
            // m odd
        {
            Result *= cos(m*phi);
            int const power = absm+2*m+1;
            Check(power%2==0);
            if ((power/2)%2) Result = -Result;
        }
        else
            // m even
        {
            Result *= sin(m*phi);
            int const power = absm+2*m+2;
            Check(power%2==0);
            if ((power/2)%2) Result = -Result;
        }
        Result *= sqrt(2.0);
    }
    else if (m>0)
    {
        if (m%2)
            // m odd
        {
            Result *= sin(m*phi);
            int const power = absm+2*m+1;
            Check(power%2==0);
            if ((power/2)%2) Result = -Result;
        }
        else
            // m even
        {
            Result *= cos(m*phi);
            int const power = absm+2*m;
            Check(power%2==0);
            if ((power/2)%2) Result = -Result;
        }
        Result *= sqrt(2.0);
    }
    return Result;
}

} // end namespace rtt_sf

//---------------------------------------------------------------------------//
//                 end of Ylm.cc
//---------------------------------------------------------------------------//
