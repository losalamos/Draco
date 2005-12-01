//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature.i.hh
 * \author Kelly Thompson
 * \date   2005 November 2005
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_quadrature_i_hh
#define rtt_quadrature_i_hh

#include <limits>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Guass-Legendre quadrature.
 *
 * Calculate abscissae and weights for Gauss-Legendre quadrature:
 * \f$\int_{0}^{\infty}{x^\alpha e^{-x} f(x) dx} = \sum_{j=0}^{N-1}w_j f(x_j)\f$
 *
 * \arg \a FieldVector A random access container on a field type.
 *
 * \param x1 Start of integration interval
 * \param x2 End of integration interval
 * \param mu The integration points.
 * \param wt The integration weights.
 * \param numGaussPoints Number of integration points in the quadrature
 * \return void but also provides mu and wt.
 */
template< typename FieldVector >
void Quadrature::calculateGaussPointsAndWeights(
    double const mu1,  // expect FieldVector::value_type to be promoted to double.
    double const mu2,
    FieldVector & mu,
    FieldVector & wt,
    unsigned const numGaussPoints )
{
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;
    using std::cos;
    using std::numeric_limits;
    
    typedef typename FieldVector::value_type Field;

    Require( numGaussPoints > 0 );
    Require( mu2 > mu1 );

    // convergence tolerance
    Field const tolerance( 100*numeric_limits< Field >::epsilon() );
    
    // size the member data vectors
    mu.resize( numGaussPoints );
    wt.resize( numGaussPoints );

    // number of Gauss points in the half range.
    // The roots are symmetric in the interval.  We only need to search for
    // half of them.
    unsigned const numHrGaussPoints( (numGaussPoints+1)/2 );

    // mid-point of integration range
    Field const mu_m( 0.5*(mu2+mu1) );
    // length of half the integration range.
    Field const mu_l( 0.5*(mu2-mu1) );
    
    // Loop over the desired roots.
    for ( size_t iroot=0; iroot<numHrGaussPoints; ++iroot)
    {
	// Approximate the i-th root.
	Field z( cos( PI * ( iroot+0.75 ) / ( numGaussPoints+0.5 )) );

        // Temporary storage;
        Field z1, pp;
        
	do // Use Newton's method to refine the value for the i-th root.  
	{
	    Field p1( 1 );
	    Field p2( 0 );

	    // This loop represents the recurrence relation needed to obtain
	    // the Legendre polynomial evaluated at z.
	    for( unsigned j=0; j<numGaussPoints; ++j )
	    {
		Field const p3 = p2;
		p2 = p1;
		p1 = (( 2*j+1) * z * p2 - j*p3 ) / (j+1);
	    }
	    
	    // p1 is now the desired Legendre polynomial evaluated at z. We
	    // next compute pp, its derivative, by a standard relation
	    // involving also p2, the polynomial of one lower order.
	    pp = numGaussPoints * (z*p1-p2)/(z*z-1.0);
	    z1 = z;
	    
	    // Newton's Method
	    z = z1 - p1/pp;

	} while( ! soft_equiv(z,z1, tolerance ) );

	// Roots will be between -1 and 1.0 and symmetric about the origin. 
        int const idxSymPart( numGaussPoints - iroot - 1 );

        // Now, scale the root to tthe desired interval and put in its
        // symmetric counterpart.
	mu[ iroot ]      = mu_m - mu_l*z;       
	mu[ idxSymPart ] = mu_m + mu_l*z;       

	// Compute the associated weight and its symmetric counterpart.
        wt[ iroot ]      = 2*mu_l / ((1-z*z)*pp*pp); 
        wt[ idxSymPart ] = wt[ iroot ];
    }	
    
    return;
}

} // end namespace rtt_quadrature

#endif // rtt_quadrature_i_hh


//---------------------------------------------------------------------------//
//                              end of Quadrature.i.hh
//---------------------------------------------------------------------------//
