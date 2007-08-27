//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Ordinate.cc
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Implementation file for the class rtt_quadrature::Ordinate.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

#include "Quadrature.hh"
#include "Ordinate.hh"
#include "GeneralQuadrature.hh"

namespace rtt_quadrature
{

using namespace std;
using rtt_dsxx::SP;


//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for OrdinateSet.
 * \param quadrature Quadrature from which to generate ordinate set
 * \param geometry Geometry of the problem.
 * \param dimension The dimension of the problem (1 or 2)
 *
 * \pre \c quadrature!=SP<Quadrature>()
 *
 * \note The quadrature object must be a level set quadrature, and it must
 * supply the number of levels.  At present all we can do is *assume* it is a
 * level set (since there is presently no way to query the object) and that
 * the number of levels equals the Sn order.
 * \todo The insertion of starting ordinates uses an algorithm that is \f$L^3\f$
 * in the number of levels \f$L\f$.  This could conceivably bite us someday
 * if computing power becomes great enough for computations with very large
 * \f$L\f$. 
 *
 * \note Notice it defaults to 2 spatial dimensions. See class boundary_source
 * for a reason why.
 */

OrdinateSet::OrdinateSet( SP<Quadrature const>       const quadrature_,
                          rtt_mesh_element::Geometry const geometry_,
                          unsigned                   const dimension_ )
    : quadrature( quadrature_ ),
      geometry(   geometry_   ),
      dimension(  dimension_  )
{
    Require( quadrature!=SP<Quadrature>()               );
    Require( quadrature->dimensionality() == 1 ||
             quadrature->dimensionality() == 2          );
    Require( dimension == 1 || dimension == 2 );

    // vector<Ordinate> ordinates;

    if( quadrature->dimensionality() == 1 ) 
        create_set_from_1d_quadrature();
    if( quadrature->dimensionality() == 2 && dimension == 2 )
        create_set_from_2d_quadrature_for_2d_mesh();
    if( quadrature->dimensionality() == 2 && dimension == 1 )
        create_set_from_2d_quadrature_for_1d_mesh();

    // OrdinateSet does not support 3D quadrature sets at this time.
    Ensure( quadrature->dimensionality() < 3 );
    Ensure( ordinates.size() > 0 );    
}
    
//---------------------------------------------------------------------------//
/*! 
 * \brief STL-compatible comparator predicate to sort ordinate vectors by xi then mu.
 * \param a
 * First comparand
 * \param b
 * Second comparand
 * \return \c true if the first comparand has a smaller xi than the second,
 * or if the xis are equal and the first comparand has a smaller mu than the
 * second; \c false otherwise.
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
    if (a.xi() < b.xi())
    {
	return true;
    }
    else if (a.xi() > b.xi())
    {
	return false;
    }
    else
    {
	return (a.mu() < b.mu());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Use the same representation as rtt_sf::galerkinYlk.
 *
 * \param l Compute spherical harmonic of degree l.
 * \param k Compute spherical harmonic of order k.
 * \param ordinate Direction cosines
 *
 * \pre \c -l<=k<=l
 * \bug depricated!!!
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

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 1D quadrature specification.
 */
void OrdinateSet::create_set_from_1d_quadrature( void )
{
    unsigned const number_of_ordinates = quadrature->getNumOrdinates();
    ordinates.resize( number_of_ordinates );
    
    for (unsigned a=0; a<number_of_ordinates; a++)
    {
        double const mu = quadrature->getMu(a);
        double const weight = quadrature->getWt(a);
        ordinates[a] = Ordinate(mu, weight);
    }
    
    if( geometry ==  rtt_mesh_element::SPHERICAL )
    {
        Insist(quadrature->dimensionality() == 1, "Quadrature dimensionality != 1");

        // insert mu=-1 starting direction 
        vector<Ordinate>::iterator a = ordinates.begin();
        a = ordinates.insert(a, Ordinate(-1.0,
                                         0.0,
                                         0.0,
                                         0.0));
    }
    else if ( geometry ==  rtt_mesh_element::CARTESIAN)
    {
        Insist(quadrature->dimensionality() == 1, "Quadrature dimensionality != 1");
    }
    else
    {
        Check(geometry == rtt_mesh_element::AXISYMMETRIC);

        Insist(false, "Axisymmetric geometry is incompatible with "
               "a 1-D quadrature set");
    }
    
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_2d_mesh( void )
{
    unsigned const number_of_ordinates = quadrature->getNumOrdinates();
    unsigned const number_of_levels = quadrature->getSnOrder();
        
    // If the geometry is axisymmetric, reserve enough room for both the
    // ordinate quadrature set and the supplemental eta==0, mu<0 starting
    // ordinates.  The latter are required to supply a starting value for
    // the ordinate differencing.
    
    ordinates.reserve(number_of_ordinates + number_of_levels);
    ordinates.resize(number_of_ordinates);
    
    // Copy the ordinates, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level
    // sets that makes it simpler to insert the supplemental ordinates
    // and set up the associated task dependencies in axisymmetric
    // geometry.
    
    if( quadrature->getEta().empty() )
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature->getMu(a);
            double const xi = quadrature->getXi(a);
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature->getWt(a);
            ordinates[a] = Ordinate(mu, eta, xi, weight);
        }
    }
    else // assume xi is empty.
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature->getMu(a);
            double const eta = quadrature->getEta(a);
            double const xi = sqrt(1-eta*eta-mu*mu);
            double const weight = quadrature->getWt(a);
            ordinates[a] = Ordinate(mu, xi, eta, weight);
        }
    }       
        
    sort( ordinates.begin(), ordinates.end(), Ordinate::SnCompare );
    
    if( geometry == rtt_mesh_element::AXISYMMETRIC )
    {
        // Define an impossible value for a direction cosine.  We use
        // this to simplify the logic of determining when we are at
        // the head of a new level set.
        
        double const SENTINEL_COSINE = 2.0;  
        
        // Insert the supplemental ordinates.  Count the levels as a sanity
        // check.
        
        unsigned check_number_of_levels = 0;
        double xi = -SENTINEL_COSINE;
        for ( vector<Ordinate>::iterator a = ordinates.begin();
              a != ordinates.end();
              ++a)
        {
            double const old_xi = xi;
            xi = a->xi();
            if (xi != old_xi)
                // We are at the start of a new level.  Insert the starting
                // ordinate.  This has eta==0 and mu determined by the
                // normalization condition.
            {
                check_number_of_levels++;
                Check(1.0-xi*xi >= 0.0);
                a = ordinates.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                                 0.0,
                                                 xi,
                                                 0.0));
            }
        }
        Check(number_of_levels==check_number_of_levels);
    }
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature
 * specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_1d_mesh( void )
{
    // Define an impossible value for a direction cosine.  We use
    // this to simplify the logic of determining when we are at
    // the head of a new level set.

    Insist(geometry == rtt_mesh_element::AXISYMMETRIC,
           "Mesh geometry != AXISYMMETRIC");
    
    double const SENTINEL_COSINE = 2.0;
    
    unsigned const number_of_ordinates = quadrature->getNumOrdinates()/2;
    unsigned const number_of_levels = quadrature->getSnOrder()/2;

    ordinates.reserve(number_of_ordinates + number_of_levels);
    ordinates.resize(number_of_ordinates);

    // Copy the ordinates, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level sets
    // that makes it simpler to insert the supplemental ordinates and set up
    // the associated task dependencies in axisymmetric geometry.
    
    unsigned check_number_of_ordinates = 0;
    for (unsigned a=0; a<2*number_of_ordinates; a++)
    {
        double const mu = quadrature->getMu(a);
        double const xi = quadrature->getXi(a);
        
        // \todo Here we check for ordinates only for \f$\xi > 0\f$ because
        // we are reducing the 2D quadrature to 1D cylindrical geometry which
        // needs quadrature ordinates only in the octant with
        // \f$\xi > 0\f$ and \f$\eta > 0\f$ \f$\mu \in [-1,1]\f$.
        // Again, logic for this should be included in Ordinate.cc.
        
        if (xi >= 0)
        {
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature->getWt(a);
            ordinates[check_number_of_ordinates] =
                Ordinate(mu, eta, xi, weight);
            ++check_number_of_ordinates;
        }
    }
    Check(number_of_ordinates==check_number_of_ordinates);
    
    sort( ordinates.begin(), ordinates.end(), Ordinate::SnCompare);
    
    // Insert the supplemental ordinates.  Count the levels as a sanity
    // check.
    
    unsigned check_number_of_levels = 0;
    double xi = -SENTINEL_COSINE;
    for( vector<Ordinate>::iterator a = ordinates.begin();
         a != ordinates.end();
         ++a )
    {
        double const old_xi = xi;
        xi = a->xi();
        if (xi != old_xi)
            // We are at the start of a new level.  Insert the starting
            // ordinate.  This has eta==0 and mu determined by the
            // normalization condition.
        {
            check_number_of_levels++;
            Check(1.0-xi*xi >= 0.0);
            a = ordinates.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                             0.0,
                                             xi,
                                             0.0));
        }
    }
    Check(number_of_levels==check_number_of_levels);
    return;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Ordinate.cc
//---------------------------------------------------------------------------//
