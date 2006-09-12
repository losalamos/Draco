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
using rtt_quadrature::Quadrature;

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for OrdinateSet.
 * \param quadrature Quadrature from which to generate ordinate set
 * \param geometry Geometry of the problem.
 * \param mesh_dimension The dimension of the mesh (1 or 2).
 *
 * \pre \c quadrature!=SP<Quadrature>()
 *
 * \note The quadrature object must be a level set quadrature, and it must
 * supply the number of levels.  At present all we can do is *assume* it is a
 * level set (since there is presently no way to query the object) and that
 * the number of levels equals the Sn order.
 * \todo The insertion of starting angles uses an algorithm that is \f$L^3\f$
 * in the number of levels \f$L\f$.  This could conceivably bite us someday
 * if computing power becomes great enough for computations with very large
 * \f$L\f$. 
 *
 * \note Notice it defaults to 2 spatial dimensions. See class boundary_source
 * for a reason why.
 */

OrdinateSet::OrdinateSet( SP<Quadrature const>       const quadrature_,
                          rtt_mesh_element::Geometry const geometry_,
                          unsigned                   const mesh_dimension_ )
    : quadrature(     quadrature_     ),
      geometry(       geometry_       ),
      mesh_dimension( mesh_dimension_ )
{
    Require( quadrature!=SP<Quadrature>()               );
    Require( quadrature->dimensionality() == 1 ||
             quadrature->dimensionality() == 2          );
    Require( mesh_dimension == 1 || mesh_dimension == 2 );

    // vector<Ordinate> angles;

    if( quadrature->dimensionality() == 1 ) 
        create_set_from_1d_quadrature();
    if( quadrature->dimensionality() == 2 && mesh_dimension == 2 )
        create_set_from_2d_quadrature_for_2d_mesh();
    if( quadrature->dimensionality() == 2 && mesh_dimension == 1 )
        create_set_from_2d_quadrature_for_1d_mesh();

    // OrdinateSet does not support 3D quadrature sets at this time.
    Ensure( quadrature->dimensionality() < 3 );
    Ensure( size() > 0 );    
}

// OrdinateSet::OrdinateSet( Quadrature                 const *quadrature_,
//                           rtt_mesh_element::Geometry const geometry_,
//                           unsigned                   const mesh_dimension_ )
//     : quadrature(     SP<Quadrature const>(quadrature_) ),
//       geometry(       geometry_       ),
//       mesh_dimension( mesh_dimension_ )
// {
//     Require( quadrature!=SP<Quadrature>()               );
//     Require( quadrature->dimensionality() == 1 ||
//              quadrature->dimensionality() == 2          );
//     Require( mesh_dimension == 1 || mesh_dimension == 2 );

//     // vector<Ordinate> angles;

//     if( quadrature->dimensionality() == 1 ) 
//         create_set_from_1d_quadrature();
//     if( quadrature->dimensionality() == 2 && mesh_dimension == 2 )
//         create_set_from_2d_quadrature_for_2d_mesh();
//     if( quadrature->dimensionality() == 2 && mesh_dimension == 1 )
//         create_set_from_2d_quadrature_for_1d_mesh();

//     // OrdinateSet does not support 3D quadrature sets at this time.
//     Ensure( quadrature->dimensionality() < 3 );
//     Ensure( size() > 0 );    
// }

    
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
 * vector< Ordinate > angles;
 * for( int i=0; i<numAngles; ++i )
 *   angles.push_back( Ordinate( spQ->getMu(i), spQ->getWt(i) ) );
 * sort(angles.begin(), angles.end(), OrdinateSet::SnCompare );
 * \endcode
 */
bool OrdinateSet::SnCompare(Ordinate const &a, Ordinate const &b)
{
    // Note that x==r==mu, z==xi
    if (a.z < b.z)
    {
	return true;
    }
    else if (a.z > b.z)
    {
	return false;
    }
    else
    {
	return (a.x < b.x);
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
double OrdinateSet::Y( unsigned const l,
                       int      const k,
                       Ordinate const &ordinate,
                       double   const sumwt )
{
    Require(abs(k)<=l);
    
    // double const theta = acos(ordinate.z);
    double const phi = atan2(ordinate.y, ordinate.x);
    return rtt_sf::galerkinYlk(l, k, ordinate.z, phi, sumwt );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 1D quadrature specification.
 */
void OrdinateSet::create_set_from_1d_quadrature( void )
{
    unsigned const number_of_angles = quadrature->getNumAngles();
    resize( number_of_angles );
    
    for (unsigned a=0; a<number_of_angles; a++)
    {
        double const mu = quadrature->getMu(a);
        double const weight = quadrature->getWt(a);
        operator[](a) = Ordinate(mu, weight);
    }
    
    if( geometry ==  rtt_mesh_element::SPHERICAL )
    {
        Insist(quadrature->dimensionality() == 1, "Quadrature dimensionality != 1");

        // insert mu=-1 starting direction 
        OrdinateSet::iterator a = begin();
        a = insert(a, Ordinate(-1.0,
                               0.0,
                               0.0,
                               0.0));
    }
    else if ( geometry ==  rtt_mesh_element::CARTESIAN)
    {
        Insist(quadrature->dimensionality() == 1, "Quadrature dimensionality != 1");
    }
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_2d_mesh( void )
{
    unsigned const number_of_angles = quadrature->getNumAngles();
    unsigned const number_of_levels = quadrature->getSnOrder();
        
    // If the geometry is axisymmetric, reserve enough room for both the
    // angle quadrature set and the supplemental eta==0, mu<0 starting
    // angles.  The latter are required to supply a starting value for
    // the angle differencing.
    
    reserve(number_of_angles + number_of_levels);
    resize(number_of_angles);
    
    // Copy the angles, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level
    // sets that makes it simpler to insert the supplemental angles
    // and set up the associated task dependencies in axisymmetric
    // geometry.
    
    if( quadrature->getEta().empty() )
    {
        for (unsigned a=0; a<number_of_angles; a++)
        {
            double const mu = quadrature->getMu(a);
            double const xi = quadrature->getXi(a);
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature->getWt(a);
            operator[](a) = Ordinate(mu, eta, xi, weight);
        }
    }
    else // assume xi is empty.
    {
        for (unsigned a=0; a<number_of_angles; a++)
        {
            double const mu = quadrature->getMu(a);
            double const eta = quadrature->getEta(a);
            double const xi = sqrt(1-eta*eta-mu*mu);
            double const weight = quadrature->getWt(a);
            operator[](a) = Ordinate(mu, xi, eta, weight);
        }
    }       
        
    sort( begin(), end(), SnCompare );
    
    if( geometry == rtt_mesh_element::AXISYMMETRIC )
    {
        // Define an impossible value for a direction cosine.  We use
        // this to simplify the logic of determining when we are at
        // the head of a new level set.
        
        double const SENTINEL_COSINE = 2.0;  
        
        // Insert the supplemental angles.  Count the levels as a sanity
        // check.
        
        unsigned check_number_of_levels = 0;
        double xi = -SENTINEL_COSINE;
        for ( iterator a = begin(); a != end(); ++a)
        {
            double const old_xi = xi;
            xi = a->z;
            if (xi != old_xi)
                // We are at the start of a new level.  Insert the starting
                // angle.  This has eta==0 and mu determined by the
                // normalization condition.
            {
                check_number_of_levels++;
                Check(1.0-xi*xi >= 0.0);
                a = insert(a, Ordinate(-sqrt(1.0-xi*xi),
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
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_1d_mesh( void )
{
    // Define an impossible value for a direction cosine.  We use
    // this to simplify the logic of determining when we are at
    // the head of a new level set.

    Insist(geometry == rtt_mesh_element::AXISYMMETRIC, "Mesh geometry != AXISYMMETRIC");
    
    double const SENTINEL_COSINE = 2.0;
    
    unsigned const number_of_angles = quadrature->getNumAngles()/2;
    unsigned const number_of_levels = quadrature->getSnOrder()/2;

    reserve(number_of_angles + number_of_levels);
    resize(number_of_angles);

    // Copy the angles, then sort -- first by xi (into level sets) and second
    // by mu.  This yields a consistent structure for the level sets that
    // makes it simpler to insert the supplemental angles and set up the
    // associated task dependencies in axisymmetric geometry.
    
    unsigned check_number_of_angles = 0;
    for (unsigned a=0; a<2*number_of_angles; a++)
    {
        double const mu = quadrature->getMu(a);
        double const xi = quadrature->getXi(a);
        
        // \todo Here we check for angles only for \f$\xi > 0\f$ because
        // we are reducing the 2D quadrature to 1D cylindrical geometry which
        // needs quadrature angles only in the octant with
        // \f$\xi > 0\f$ and \f$\eta > 0\f$ \f$\mu \in [-1,1]\f$.
        // Again, logic for this should be included in Ordinate.cc.
        
        if (xi >= 0)
        {
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature->getWt(a);
            operator[](check_number_of_angles) = Ordinate(mu, eta, xi, weight);
            ++check_number_of_angles;
        }
    }
    Check(number_of_angles==check_number_of_angles);
    
    sort( begin(), end(), SnCompare);
    
    // Insert the supplemental angles.  Count the levels as a sanity
    // check.
    
    unsigned check_number_of_levels = 0;
    double xi = -SENTINEL_COSINE;
    for( iterator a = begin(); a != end(); ++a )
    {
        double const old_xi = xi;
        xi = a->z;
        if (xi != old_xi)
            // We are at the start of a new level.  Insert the starting
            // angle.  This has eta==0 and mu determined by the
            // normalization condition.
        {
            check_number_of_levels++;
            Check(1.0-xi*xi >= 0.0);
            a = insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                   0.0,
                                   xi,
                                   0.0));
        }
    }
    Check(number_of_levels==check_number_of_levels);
    return;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Create a new quadrature set that has the starting directions.
 */
// rtt_dsxx::SP< const Quadrature > OrdinateSet::getQuadrature( void ) const
// {
//     // We only need an augmented quadrature set if we are are using
//     // axisymmetric geometry.
//     if( geometry != rtt_mesh_element::AXISYMMETRIC )
//         return quadrature;

//     size_t len( size() );
//     std::vector<double> mu( len );
//     std::vector<double> eta(len );
//     std::vector<double> xi( len );
//     std::vector<double> wt( len );
//     for( size_t i=0; i<len; ++i )
//     {
//         Ordinate omega=this->operator[](i);
//         mu[i]  = omega.mu();
//         eta[i] = omega.eta();
//         xi[i]  = omega.xi();
//         wt[i]  = omega.wt();
//     }
//     std::string qname("Augmented axisymmetric ");
//     qname += quadrature->name();
    
//     rtt_dsxx::SP< const Quadrature > spQuad(
//         new GeneralQuadrature( quadrature->getSnOrder(),
//                                quadrature->getNorm(),
//                                mu,
//                                eta,
//                                xi,
//                                wt,
//                                quadrature->getSnOrder(),
//                                quadrature->dimensionality(),
//                                qname ) );
//     return spQuad;
// }


} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Ordinate.cc
//---------------------------------------------------------------------------//
