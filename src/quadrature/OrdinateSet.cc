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
// $Id: Ordinate.cc 6499 2012-03-15 20:19:33Z kgbudge $
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

#include "Quadrature.hh"
#include "GeneralQuadrature.hh"
#include "OrdinateSet.hh"
#include "QuadServices_GQ.hh"
#include "QuadServices_SN.hh"

namespace rtt_quadrature
{
using namespace std;
using namespace rtt_dsxx;

std::vector<Ordinate> OrdinateSet::set_ordinates(rtt_mesh_element::Geometry const geometry)
{
    std::vector<Ordinate> Result;

    if( quadrature_->dimensionality() == 1 )
    {
        if (dimension_ == 1)
        {
            Result = create_set_from_1d_quadrature(geometry);
        }
        else if (dimension_ == 2)
        {
            throw invalid_argument("cannot construct 2D ordinate set "
                                   "from 1D quadrature");
        }
        else
        {
            Check(dimension_ == 3);
            throw invalid_argument("cannot construct 3D ordinate set "
                                   "from 1D quadrature");
        }
    }
    else if( quadrature_->dimensionality() == 2)
    {
        if (dimension_ == 1 )
        {
            Result = create_set_from_2d_quadrature_for_1d_mesh(geometry);
        }
        else if (dimension_ == 2 )
        {
            Result = create_set_from_2d_quadrature_for_2d_mesh(geometry);
        }
        else
        {
            Check(dimension_ == 3 );
            Result = create_set_from_2d_quadrature_for_3d_mesh();
        }
    }
    else
    {
        Check(quadrature_->dimensionality() == 3);
        if (dimension_ == 1)
        {
            throw invalid_argument("sorry, not implemented: 1D ordinate set "
                                   "from 3D quadrature");
        }
        else if (dimension_ == 2)
        {
            throw invalid_argument("sorry, not implemented: 2D ordinate set "
                                   "from 3D quadrature");
        }
        else
        {
            Check(dimension_ == 3);
            Result = create_set_from_3d_quadrature_for_3d_mesh();
        }
    }

    // Ensure that the norm is the same as that of the quadrature.
    double norm = 0;
    unsigned const N = Result.size();
    for (unsigned a=0; a<N; ++a)
    {
        norm += Result[a].wt();
    }
    norm_ = quadrature_->getNorm();
    norm = norm_/norm;
    if (norm != 1.0)
    {
        for (unsigned a=0; a<N; ++a)
        {
            Result[a] = Ordinate(Result[a].mu(),
                                 Result[a].eta(),
                                 Result[a].xi(),
                                 Result[a].wt()*norm);
        }
    }

    std::sort( Result.begin(), Result.end(), comparator_); 

    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for OrdinateSet.
 *
 * \param quadrature Quadrature from which to generate ordinate set
 * \param geometry Geometry of the problem.
 * \param dimension The dimension of the problem
 * \param extra_starting_directions Should be used only with curvilinear
 * geometry. Specifies that the R-reflection of each starting direction should
 * be included in the ordinate set.
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
                          unsigned                   const expansion_order,
                          Quadrature::QIM            const qim,
                          bool                       const extra_starting_directions,
                          comparator_t               const comparator) 
    : quadrature_( quadrature ),
      geometry_(   geometry   ),
      dimension_(  dimension  ),
      expansion_order_(  expansion_order  ),
      extra_starting_directions_( extra_starting_directions ),
      comparator_( comparator ),
      ordinates_(set_ordinates(geometry))
{
    Require( quadrature!=SP<Quadrature>()               );
    Require( quadrature->dimensionality() == 1 ||
             quadrature->dimensionality() == 2 ||
             quadrature->dimensionality() == 3          );
    Require( dimension == 1 || dimension == 2 || dimension == 3);
    Require( (geometry==rtt_mesh_element::SPHERICAL && dimension == 1 )  ||
             (geometry==rtt_mesh_element::AXISYMMETRIC && dimension <3 ) ||
             geometry==rtt_mesh_element::CARTESIAN);

    // This version sets the interpolation model based on a construction parameter.

    if (qim == Quadrature::SN)
    {
        qs_ = SP<QuadServices>(new QuadServices_SN(ordinates_, norm_, dimension_, expansion_order_, geometry));
    }
    else if (qim == Quadrature::GQ)
    { 
        qs_ = SP<QuadServices>(new QuadServices_GQ(ordinates_, norm_, dimension_, expansion_order_, geometry));
    }

    Ensure( check_class_invariants() );
    Ensure( getNumOrdinates() > 0 );
    Ensure( getQuadrature() == quadrature );
    Ensure( getGeometry() == geometry );
    Ensure( getDimension() == dimension );
}

OrdinateSet::OrdinateSet( SP<Quadrature const>       const quadrature,
                          rtt_mesh_element::Geometry const geometry,
                          unsigned                   const dimension,
                          unsigned                   const expansion_order,
                          bool                       const extra_starting_directions,
                          comparator_t               const comparator) 
    : quadrature_( quadrature ),
      geometry_(   geometry   ),
      dimension_(  dimension  ),
      expansion_order_(  expansion_order  ),
      extra_starting_directions_( extra_starting_directions ),
      comparator_( comparator ),
      ordinates_(set_ordinates(geometry))
{
    Require( quadrature!=SP<Quadrature>()               );
    Require( quadrature->dimensionality() == 1 ||
             quadrature->dimensionality() == 2 ||
             quadrature->dimensionality() == 3          );
    Require( dimension == 1 || dimension == 2 || dimension == 3);
    Require( (geometry==rtt_mesh_element::SPHERICAL && dimension == 1 )  ||
             (geometry==rtt_mesh_element::AXISYMMETRIC && dimension <3 ) ||
             geometry==rtt_mesh_element::CARTESIAN);

    // This version sets the interpolation model based on the quadrature object.

    if (quadrature_->interpolation_model() == Quadrature::SN)
    {
        qs_ = SP<QuadServices>(new QuadServices_SN(ordinates_, norm_, dimension_, expansion_order_, geometry));
    }
    else if (quadrature_->interpolation_model() == Quadrature::GQ)
    { 
        qs_ = SP<QuadServices>(new QuadServices_GQ(ordinates_, norm_, dimension_, expansion_order_, geometry));
    }


//-----------------------------------------------

    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << "ORDINATES" << endl << endl;
    cout << "   m  \t     mu        \t     eta       \t      xi       \t      Phi      \t      wt      " << endl;
    cout << "  --- \t-------------- \t-------------- \t-------------- \t-------------- \t--------------" << endl;
    double sum_wt = 0.0;
    for ( size_t m = 0; m < ordinates_.size(); ++m)
    {
        double const phi((dimension_== 1 && geometry_ == rtt_mesh_element::CARTESIAN)
                         ? 0
                         : qs_->compute_azimuthalAngle(ordinates_[m].mu(), ordinates_[m].eta(), ordinates_[m].xi()));
        cout << "   "
	     << m << "\t"
	     << setprecision(10) << setw(14) << ordinates_[m].mu()  << "\t"
	     << setprecision(10) << setw(14) << ordinates_[m].eta() << "\t"
	     << setprecision(10) << setw(14) << ordinates_[m].xi()  << "\t"
	     << setprecision(4)  << setw(8)  << phi << "(" <<  setw(6)  << phi*360/2/rtt_units::PI << ")" << "\t" 
	     << setprecision(10) << setw(14) << ordinates_[m].wt()  << endl;
	sum_wt += ordinates_[m].wt();
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;

//-----------------------------------------------

    Ensure( check_class_invariants() );
    Ensure( getNumOrdinates() > 0 );
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
/*!\brief Helper for creating an OrdinateSet from a 1D quadrature
 * specification.
 */
std::vector<Ordinate> 
OrdinateSet::create_set_from_1d_quadrature(rtt_mesh_element::Geometry const geometry)
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    std::vector<Ordinate> Result( number_of_ordinates );
    
    for (unsigned a=0; a<number_of_ordinates; a++)
    {
        double const mu = quadrature_->getMu(a);
        double const weight = quadrature_->getWt(a);
        Result[a] = Ordinate(mu, weight);
    }
    
    if( geometry ==  rtt_mesh_element::SPHERICAL )
    {
        Insist(quadrature_->dimensionality() == 1,
               "Quadrature dimensionality != 1");

        // insert mu=-1 starting direction 
        vector<Ordinate>::iterator a = Result.begin();
        a = Result.insert(a, Ordinate(-1.0,
                                          0.0,
                                          0.0,
                                          0.0));
        
        // insert mu=1 starting direction 
        if (extra_starting_directions_)
            Result.push_back(Ordinate(1.0,
                                          0.0,
                                          0.0,
                                          0.0));
    }
    else if ( geometry ==  rtt_mesh_element::CARTESIAN)
    {
        Insist(quadrature_->dimensionality() == 1, "Quadrature dimensionality != 1.");
    }
    else
    {
        Check(geometry == rtt_mesh_element::AXISYMMETRIC);

        Insist(false, "Axisymmetric geometry is incompatible with 1-D quadrature sets.");
    }
    
    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
std::vector<Ordinate> 
OrdinateSet::create_set_from_2d_quadrature_for_2d_mesh(rtt_mesh_element::Geometry const geometry)
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    unsigned const number_of_levels = quadrature_->getSnOrder();
        
    // If the geometry is axisymmetric, reserve enough room for both the
    // ordinate quadrature set and the supplemental eta==0, mu<0 starting
    // ordinates.  The latter are required to supply a starting value for
    // the ordinate differencing.
    
    std::vector<Ordinate> Result;
    Result.reserve(number_of_ordinates + number_of_levels);
    Result.resize(number_of_ordinates);
    
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
            Result[a] = Ordinate(mu, eta, xi, weight);
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
            Result[a] = Ordinate(mu, xi, eta, weight);
        }
    }       

    sort( Result.begin(), Result.end(), comparator_);

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

        for ( vector<Ordinate>::iterator a = Result.begin(); a != Result.end(); ++a)
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
                a = Result.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                              0.0,
                                              xi,
                                              0.0));
                
                // insert mu > 0
                if (extra_starting_directions_)
                    if (a != Result.begin())
                        a = Result.insert(a, Ordinate(sqrt(1.0-old_xi*old_xi),
                                                      0.0,
                                                      old_xi,
                                                      0.0));
            }
        }

        // insert mu > 0 on the final level
        if (extra_starting_directions_)
            Result.push_back(Ordinate(sqrt(1.0-xi*xi),
                                      0.0,
                                      xi,
                                      0.0));
        
        Check(number_of_levels==check_number_of_levels);
    }

    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 2D quadrature specification.
 */
std::vector<Ordinate> 
OrdinateSet::create_set_from_2d_quadrature_for_1d_mesh(rtt_mesh_element::Geometry const geometry)
{
    // Define an impossible value for a direction cosine.  We use
    // this to simplify the logic of determining when we are at
    // the head of a new level set.

    Insist(geometry == rtt_mesh_element::AXISYMMETRIC, "Mesh geometry is not AXISYMMETRIC");
    
    double const SENTINEL_COSINE = 2.0;
    
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates()/2;
    unsigned const number_of_levels = quadrature_->getSnOrder()/2;

    std::vector<Ordinate> Result;
    Result.reserve(number_of_ordinates + number_of_levels);
    Result.resize(number_of_ordinates);

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
            Result[check_number_of_ordinates] = Ordinate(mu, eta, xi, weight);
            ++check_number_of_ordinates;
        }
    }
    Check(number_of_ordinates==check_number_of_ordinates);
    
    sort( Result.begin(), Result.end(), Ordinate::SnCompare);
    
    // Insert the supplemental ordinates.  Count the levels as a sanity
    // check.
    
    unsigned check_number_of_levels = 0;
    double xi = -SENTINEL_COSINE;
    for( vector<Ordinate>::iterator a = Result.begin(); a != Result.end(); ++a )
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
            a = Result.insert(a, Ordinate(-sqrt(1.0-xi*xi),
                                             0.0,
                                             xi,
                                             0.0));
            // insert mu > 0
            if (extra_starting_directions_)
                if (a != Result.begin())
                    a = Result.insert(a, Ordinate(sqrt(1.0-old_xi*old_xi),
                                                      0.0,
                                                      old_xi,
                                                      0.0));
        }
    }

    // insert mu > 0 on the final level
    if (extra_starting_directions_)
        Result.push_back(Ordinate(sqrt(1.0-xi*xi),
                                      0.0,
                                      xi,
                                      0.0));
    
    Check(number_of_levels==check_number_of_levels);

    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 3D quadrature specification.
 */
std::vector<Ordinate> 
OrdinateSet::create_set_from_2d_quadrature_for_3d_mesh()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();
    std::vector<Ordinate> Result(2*number_of_ordinates);
    
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
            Result[a] = Ordinate(mu, eta, xi, weight);
            Result[a+number_of_ordinates] = Ordinate(mu, -eta, xi, weight);
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
            Result[a] = Ordinate(mu, xi, eta, weight);
            Result[a+number_of_ordinates] = Ordinate(mu, -xi, eta, weight);
        }
    }       
        
    sort( Result.begin(), Result.end(), Ordinate::SnComparePARTISN3);

    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper for creating an OrdinateSet from a 3D quadrature specification.
 */
std::vector<Ordinate> 
OrdinateSet::create_set_from_3d_quadrature_for_3d_mesh()
{
    unsigned const number_of_ordinates = quadrature_->getNumOrdinates();

    std::vector<Ordinate> Result;
    Result.reserve(number_of_ordinates);
    Result.resize(number_of_ordinates);
    
    // Copy the ordinates

    for (unsigned a=0; a<number_of_ordinates; a++)
    {
        double const mu = quadrature_->getMu(a);
        double const xi = quadrature_->getXi(a);
        double const eta = quadrature_->getEta(a);
        double const weight = quadrature_->getWt(a);
        Result[a] = Ordinate(mu, eta, xi, weight);
    }

    sort( Result.begin(), Result.end(), Ordinate::SnComparePARTISN3 );

    comparator_ = Ordinate::SnComparePARTISN3;
    
    return Result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of OrdinateSet.cc
//---------------------------------------------------------------------------//
