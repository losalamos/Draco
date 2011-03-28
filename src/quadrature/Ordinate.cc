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
#include "Ordinate.hh"
#include "GeneralQuadrature.hh"

namespace rtt_quadrature
{
using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for OrdinateSet.
 * \param quadrature Quadrature from which to generate ordinate set
 * \param geometry Geometry of the problem.
 * \param dimension The dimension of the problem
 *
 * \note The quadrature object must be a level set quadrature, and it must
 * supply the number of levels.  At present all we can do is *assume* it is a
 * level set (since there is presently no way to query the object) and that
 * the number of levels equals the Sn order.
 *
 * \todo The insertion of starting ordinates uses an algorithm that is
 * \f$L^3\f$ in the number of levels \f$L\f$.  This could conceivably bite us
 * someday if computing power becomes great enough for computations with very
 * large \f$L\f$.
 */

OrdinateSet::OrdinateSet( SP<Quadrature const>       const quadrature,
                          rtt_mesh_element::Geometry const geometry,
                          unsigned                   const dimension,
                          bool                       const extra_starting_directions)
    : ordinates_(),
      quadrature_( quadrature ),
      geometry_(   geometry   ),
      dimension_(  dimension  ),
      norm_(0.0),
      comparator_(Ordinate::SnCompare), // default comparator
      extra_starting_directions_(extra_starting_directions) 
{
    Require( quadrature!=SP<Quadrature>()               );
    Require( quadrature->dimensionality() == 1 ||
             quadrature->dimensionality() == 2 ||
             quadrature->dimensionality() == 3          );
    Require( dimension == 1 || dimension == 2 || dimension == 3);
    Require( (geometry==rtt_mesh_element::SPHERICAL && dimension == 1 )  ||
             (geometry==rtt_mesh_element::AXISYMMETRIC && dimension <3 ) ||
             geometry==rtt_mesh_element::CARTESIAN);

    // vector<Ordinate> ordinates;

    if( quadrature->dimensionality() == 1 )
    {
        if (dimension == 1)
        {
            create_set_from_1d_quadrature();
        }
        else if (dimension == 2)
        {
            throw invalid_argument("cannot construct 2D ordinate set "
                                   "from 1D quadrature");
        }
        else
        {
            Check(dimension==3);
            throw invalid_argument("cannot construct 3D ordinate set "
                                   "from 1D quadrature");
        }
    }
    else if( quadrature->dimensionality() == 2)
    {
        if (dimension == 1 )
        {
            create_set_from_2d_quadrature_for_1d_mesh();
        }
        else if (dimension == 2 )
        {
            create_set_from_2d_quadrature_for_2d_mesh();
        }
        else
        {
            Check(dimension == 3 );
            create_set_from_2d_quadrature_for_3d_mesh();
        }
    }
    else
    {
        Check(quadrature->dimensionality() == 3);
        if (dimension == 1)
        {
            throw invalid_argument("sorry, not implemented: 1D ordinate set "
                                   "from 3D quadrature");
        }
        else if (dimension == 2)
        {
            throw invalid_argument("sorry, not implemented: 2D ordinate set "
                                   "from 3D quadrature");
        }
        else
        {
            Check(dimension == 3);
            create_set_from_3d_quadrature_for_3d_mesh();
        }
    }

    // Ensure that the norm is the same as that of the quadrature.
    double norm = 0;
    unsigned const N = ordinates_.size();
    for (unsigned a=0; a<N; ++a)
    {
        norm += ordinates_[a].wt();
    }
    norm_ = quadrature->getNorm();
    norm = norm_/norm;
    if (norm != 1.0)
    {
        for (unsigned a=0; a<N; ++a)
        {
            ordinates_[a] = Ordinate(ordinates_[a].mu(),
                                     ordinates_[a].eta(),
                                     ordinates_[a].xi(),
                                     ordinates_[a].wt()*norm);
        }
    }

    Ensure( check_class_invariants() );
    Ensure( getOrdinates().size() > 0 );
    Ensure( getQuadrature() == quadrature );
    Ensure( getGeometry() == geometry );
    Ensure( getDimension() == dimension );
}

//---------------------------------------------------------------------------//
bool OrdinateSet::check_class_invariants() const
{
    return
        quadrature_ != SP<Quadrature>() &&
        ordinates_.size()>0;
}
    
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

//---------------------------------------------------------------------------//
/*!\brief Helper for creating an OrdinateSet from a 1D quadrature
 * specification.
 */
void OrdinateSet::create_set_from_1d_quadrature()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    ordinates_.resize( number_of_ordinates );
    
    for (unsigned a=0; a<number_of_ordinates; a++)
    {
        double const mu = quadrature_->getMu(a);
        double const weight = quadrature_->getWt(a);
        ordinates_[a] = Ordinate(mu, weight);
    }
    
    if( geometry_ ==  rtt_mesh_element::SPHERICAL )
    {
        Insist(quadrature_->dimensionality() == 1,
               "Quadrature dimensionality != 1");

        // insert mu=-1 starting direction 
        vector<Ordinate>::iterator a = ordinates_.begin();
        a = ordinates_.insert(a, Ordinate(-1.0,
                                          0.0,
                                          0.0,
                                          0.0));
        
        // insert mu=1 starting direction 
        if (extra_starting_directions_)
            ordinates_.push_back(Ordinate(1.0,
                                          0.0,
                                          0.0,
                                          0.0));
    }
    else if ( geometry_ ==  rtt_mesh_element::CARTESIAN)
    {
        Insist(quadrature_->dimensionality() == 1,
               "Quadrature dimensionality != 1");
    }
    else
    {
        Check(geometry_ == rtt_mesh_element::AXISYMMETRIC);

        Insist(false,
               "Axisymmetric geometry is incompatible with "
               "a 1-D quadrature set");
    }
    
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_2d_mesh()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    unsigned const number_of_levels = quadrature_->getSnOrder();
        
    // If the geometry is axisymmetric, reserve enough room for both the
    // ordinate quadrature set and the supplemental eta==0, mu<0 starting
    // ordinates.  The latter are required to supply a starting value for
    // the ordinate differencing.
    
    ordinates_.reserve(number_of_ordinates + number_of_levels);
    ordinates_.resize(number_of_ordinates);
    
    // Copy the ordinates, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level
    // sets that makes it simpler to insert the supplemental ordinates
    // and set up the associated task dependencies in axisymmetric
    // geometry.
    
    if( quadrature_->getEta().empty() )
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature_->getMu(a);
            double const xi = quadrature_->getXi(a);
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature_->getWt(a);
            ordinates_[a] = Ordinate(mu, eta, xi, weight);
        }
    }
    else // assume xi is empty.
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature_->getMu(a);
            double const eta = quadrature_->getEta(a);
            double const xi = sqrt(1-eta*eta-mu*mu);
            double const weight = quadrature_->getWt(a);
            ordinates_[a] = Ordinate(mu, xi, eta, weight);
        }
    }       

    if (geometry_ == rtt_mesh_element::CARTESIAN)
    {
        sort( ordinates_.begin(),
              ordinates_.end(),
              Ordinate::SnComparePARTISN3 );

        comparator_ = Ordinate::SnComparePARTISN3;
    }
    else if( geometry_ == rtt_mesh_element::AXISYMMETRIC )
    {
        sort( ordinates_.begin(), ordinates_.end(), Ordinate::SnCompare );

        // Define an impossible value for a direction cosine.  We use
        // this to simplify the logic of determining when we are at
        // the head of a new level set.
        
        double const SENTINEL_COSINE = 2.0;  
        
        // Insert the supplemental ordinates.  Count the levels as a sanity
        // check.
        
        unsigned check_number_of_levels = 0;
        double xi = -SENTINEL_COSINE;

        for ( vector<Ordinate>::iterator a = ordinates_.begin(); a != ordinates_.end(); ++a)
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

                // insert mu < 0
                a = ordinates_.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                                  0.0,
                                                  xi,
                                                  0.0));

                // insert mu > 0
                if (extra_starting_directions_)
                    if (a != ordinates_.begin())
                        a = ordinates_.insert(a, Ordinate(sqrt(1.0-old_xi*old_xi),
                                                          0.0,
                                                          old_xi,
                                                          0.0));
            }
        }

        // insert mu > 0 on the final level
        if (extra_starting_directions_)
            ordinates_.push_back(Ordinate(sqrt(1.0-xi*xi),
                                          0.0,
                                          xi,
                                          0.0));

        Check(number_of_levels==check_number_of_levels);
    }
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature
 * specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_1d_mesh()
{
    // Define an impossible value for a direction cosine.  We use
    // this to simplify the logic of determining when we are at
    // the head of a new level set.

    Insist(geometry_ == rtt_mesh_element::AXISYMMETRIC,
           "Mesh geometry != AXISYMMETRIC");
    
    double const SENTINEL_COSINE = 2.0;
    
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates()/2;
    unsigned const number_of_levels = quadrature_->getSnOrder()/2;

    ordinates_.reserve(number_of_ordinates + number_of_levels);
    ordinates_.resize(number_of_ordinates);

    // Copy the ordinates, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level sets
    // that makes it simpler to insert the supplemental ordinates and set up
    // the associated task dependencies in axisymmetric geometry.
    
    unsigned check_number_of_ordinates = 0;
    for (unsigned a=0; a<2*number_of_ordinates; a++)
    {
        double const mu = quadrature_->getMu(a);
        double const xi = quadrature_->getXi(a);
        
        // \todo Here we check for ordinates only for \f$\xi > 0\f$ because
        // we are reducing the 2D quadrature to 1D cylindrical geometry which
        // needs quadrature ordinates only in the octant with
        // \f$\xi > 0\f$ and \f$\eta > 0\f$ \f$\mu \in [-1,1]\f$.
        // Again, logic for this should be included in Ordinate.cc.
        
        if (xi >= 0)
        {
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature_->getWt(a);
            ordinates_[check_number_of_ordinates] = Ordinate(mu, eta, xi, weight);
            ++check_number_of_ordinates;
        }
    }
    Check(number_of_ordinates==check_number_of_ordinates);
    
    sort( ordinates_.begin(), ordinates_.end(), Ordinate::SnCompare);
    
    // Insert the supplemental ordinates.  Count the levels as a sanity
    // check.
    
    unsigned check_number_of_levels = 0;
    double xi = -SENTINEL_COSINE;
    for( vector<Ordinate>::iterator a = ordinates_.begin(); a != ordinates_.end(); ++a )
    {
        double const old_xi = xi;
        xi = a->xi();
        if (xi != old_xi)
            // We are at the start of a new level.  Insert the starting
            // ordinate with zero weight.
            // This has eta==0 and mu determined by the normalization condition.
        {
            check_number_of_levels++;
            Check(1.0-xi*xi >= 0.0);

            // insert mu < 0
            a = ordinates_.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                             0.0,
                                             xi,
                                             0.0));
            // insert mu > 0
            if (extra_starting_directions_)
                if (a != ordinates_.begin())
                    a = ordinates_.insert(a, Ordinate(sqrt(1.0-old_xi*old_xi),
                                                      0.0,
                                                      old_xi,
                                                      0.0));
        }
    }

    // insert mu > 0 on the final level
    if (extra_starting_directions_)
        ordinates_.push_back(Ordinate(sqrt(1.0-xi*xi),
                                      0.0,
                                      xi,
                                      0.0));
    
    Check(number_of_levels==check_number_of_levels);
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 3D quadrature
 * specification.
 */
void OrdinateSet::create_set_from_2d_quadrature_for_3d_mesh()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    ordinates_.resize(2*number_of_ordinates);
    
    // Copy the ordinates, then sort -- first by xi (into level sets) and
    // second by mu.  This yields a consistent structure for the level
    // sets that makes it simpler to insert the supplemental ordinates
    // and set up the associated task dependencies in axisymmetric
    // geometry.
    
    if( quadrature_->getEta().empty() )
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature_->getMu(a);
            double const xi = quadrature_->getXi(a);
            double const eta = sqrt(1-xi*xi-mu*mu);
            double const weight = quadrature_->getWt(a);
            ordinates_[a] = Ordinate(mu, eta, xi, weight);
            ordinates_[a+number_of_ordinates] = Ordinate(mu, -eta, xi, weight);
        }
    }
    else // assume xi is empty.
    {
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = quadrature_->getMu(a);
            double const eta = quadrature_->getEta(a);
            double const xi = sqrt(1-eta*eta-mu*mu);
            double const weight = quadrature_->getWt(a);
            ordinates_[a] = Ordinate(mu, xi, eta, weight);
            ordinates_[a+number_of_ordinates] = Ordinate(mu, -xi, eta, weight);
        }
    }       
        
    sort( ordinates_.begin(), ordinates_.end(), Ordinate::SnCompare );

    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 3D quadrature
 * specification.
 */
void OrdinateSet::create_set_from_3d_quadrature_for_3d_mesh()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    
    ordinates_.reserve(number_of_ordinates);
    ordinates_.resize(number_of_ordinates);
    
    // Copy the ordinates

    for (unsigned a=0; a<number_of_ordinates; a++)
    {
        double const mu = quadrature_->getMu(a);
        double const xi = quadrature_->getXi(a);
        double const eta = quadrature_->getEta(a);
        double const weight = quadrature_->getWt(a);
        ordinates_[a] = Ordinate(mu, eta, xi, weight);
    }

    sort( ordinates_.begin(), ordinates_.end(), Ordinate::SnComparePARTISN3 );

    comparator_ = Ordinate::SnComparePARTISN3;
    
    return;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Ordinate.cc
//---------------------------------------------------------------------------//
