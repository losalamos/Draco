//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate.cc
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Implementation file for the class rtt_quadrature::Ordinate.
 * \note   Copyright ©  2006-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <iostream>

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

#include "Quadrature.hh"
#include "GeneralQuadrature.hh"

#include "Ordinate.hh"

namespace rtt_quadrature
{
using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
/*! 
 *
 * \param a
 * First comparand
 * \param b
 * Second comparand
 * \return \c true if the first comparand has a smaller xi than the second,
 * or if the xis are equal and the first comparand has a smaller mu than the
 * second, of if the xis and mus are equal and the first comparand has a
 * smaller eta than the second; \c false otherwise.
 *
 * Typical usage:
 * \code
 * vector< Ordinate > ordinates;
 * for( int i=0; i<numOrdinates; ++i )
 *   ordinates.push_back( Ordinate( spQ->getMu(i), spQ->getWt(i) ) );
 * sort(ordinates.begin(), ordinates.end(), Ordinate::SnCompare );
 * \endcode
 */
bool Ordinate::SnCompare(Ordinate const &a, Ordinate const &b)
{
    // Note that x==r==mu, z==xi

    if (soft_equiv(a.xi(), b.xi()) && soft_equiv(a.mu(), b.mu()) && soft_equiv(a.eta(), b.eta()) )
    {
        return false;
    }
    else if (a.xi() < b.xi()) 
    {
	return true;
    }
    else if (a.xi() > b.xi()) 
    {
	return false;
    }
    else if (a.mu() < b.mu())
    {
	return true;
    }
    else if (a.mu() > b.mu())
    {
        return false;
    }
    else
    {
        return a.eta() < b.eta();
    }
}
  
//---------------------------------------------------------------------------//
/*! 
 *
 * \param a
 * First comparand
 * \param b
 * Second comparand
 * \return \c true if the first comparand is lower in PARTISN-3D order than
 * the second comparand. PARTISN order is by octant first, further ordered by
 * sign of xi, then sign of eta, then sign of mu, and then by absolute value
 * of xi, then eta, then mu within each octant.
 *
 * Typical usage:
 * \code
 * vector< Ordinate > ordinates;
 * for( int i=0; i<numOrdinates; ++i )
 *   ordinates.push_back( Ordinate( spQ->getMu(i), spQ->getWt(i) ) );
 * sort(ordinates.begin(), ordinates.end(), Ordinate::SnComparePARTISN3 );
 * \endcode
 */
bool Ordinate::SnComparePARTISN3(Ordinate const &a, Ordinate const &b)
{
    // Note that x==r==mu, z==xi
    //if (soft_equiv(a.xi(), b.xi()) && soft_equiv(a.mu(), b.mu()) && soft_equiv(a.eta(), b.eta()) )
    //{
    //    return false;
    //} else
    if (a.xi()<0 && b.xi()>0)
    {
        return true;
    }
    else if (a.xi()>0 && b.xi()<0)
    {
        return false;
    }
    else if (a.eta()<0 && b.eta()>0)
    {
        return true;
    }
    else if (a.eta()>0 && b.eta()<0)
    {
        return false;
    }
    else if (a.mu()<0 && b.mu()>0)
    {
        return true;
    }
    else if (a.mu()>0 && b.mu()<0)
    {
        return false;
    }
    else if (!soft_equiv(fabs(a.xi()), fabs(b.xi()), 1.0e-14))
    {
        return (fabs(a.xi()) < fabs(b.xi()));
    }
    else if (!soft_equiv(fabs(a.eta()), fabs(b.eta()), 1.0e-14))
    {
        return (fabs(a.eta()) < fabs(b.eta()));
    }
    else
    {
	return (!soft_equiv(fabs(a.mu()), fabs(b.mu()), 1.0e-14) &&
                fabs(a.mu()) < fabs(b.mu()));
    }
}

//---------------------------------------------------------------------------//
/*!
 * This function uses the same representation as rtt_sf::galerkinYlk.
 *
 * \param l Compute spherical harmonic of degree l.
 * \param k Compute spherical harmonic of order k.
 * \param ordinate Direction cosines
 */
double Ordinate::Y( unsigned const l,
                    int      const k,
                    Ordinate const &ordinate,
                    double   const sumwt )
{
    Require(static_cast<unsigned>(abs(k))<=l);
    
    // double const theta = acos(ordinate.xi());
    double const phi = atan2(ordinate.eta(), ordinate.mu());
    return rtt_sf::galerkinYlk(l, k, ordinate.xi(), phi, sumwt );
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Ordinate.cc
//---------------------------------------------------------------------------//
