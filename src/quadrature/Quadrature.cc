//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 15:38:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// Revision History
// -----------------
// 1.5) Moved "using" statements inside of namespace block.
//      Moved declaration of EPS inside of namespace block.
//      Fixed ensure statements to take into account the possible
//         normalization of the quadrature weights.
//      General code cleanup.

#include "Quadrature.hh"

namespace rtt_quadrature
{

// BASE CLASS QUADRATURE MEMBER FUNCTION DEFINITIONS

// The member functions included in the base class can be called for all
// quadrature sets.
//------------------------------------------------------------------------------

// Functions that return quadrature directions and weights.

// vector<double> Quadrature::getOmega( const int m ) const {
//     vector<double> omega(3);
//     Require( m >= 0 && m < getNumAngles() );
//     omega[0] = mu[m];
//     if ( dimensionality() >= 2 )
// 	omega[1] = eta[m];
//     else
// 	omega[1] = 0.0;  // no eta or xi values in 1D!
//     if ( dimensionality() >= 3 )
// 	omega[2] = xi[m];
//     else
// 	omega[2] = 0.0;  // no eta or xi values in 1D!
//     return omega;
// }

// Functions that test the validity of the quadrature set.

//---------------------------------------------------------------------------//


/*!
 * \brief Integrates dOmega over the unit sphere. (The sum of quadrature weights.)
 */
double Quadrature::iDomega() const {
    double integral = 0.0;
    for ( size_t i = 0; i < getNumAngles(); ++i )
	integral += wt[i];
    return integral;
}

//---------------------------------------------------------------------------//

/*!
 * \brief Integrates the vector Omega over the unit sphere. 
 *
 * The solution to this integral is a vector with length equal to the
 * number of dimensions of the quadrature set.  The solution should be
 * the zero vector. 
 *
 * The integral is actually calculated as a quadrature sum over all
 * directions of the quantity: 
 *
 *     Omega(m) * wt(m)
 *
 * Omega is a vector quantity.
 */
vector<double> Quadrature::iOmegaDomega() const {
    size_t ndims = dimensionality();
    vector<double> integral(ndims);
    // initialize the sum to zero.
    for ( size_t j = 0; j < ndims; ++j ) 
	integral[j] = 0.0;
    switch( ndims ) {
    case 3:
	for ( size_t i = 0; i < getNumAngles(); ++i )
	    integral[2] += wt[i]*xi[i];
	//lint -fallthrough
    case 2:
	for ( size_t i = 0; i < getNumAngles(); ++i )
	    integral[1] += wt[i]*eta[i];
	//lint -fallthrough
    case 1:
	for ( size_t i = 0; i < getNumAngles(); ++i )
	    integral[0] += wt[i]*mu[i];
	break;
    default:
	Insist(false,"Number of spatial dimensions must be 1, 2 or 3.");
    }
    return integral;
}

//---------------------------------------------------------------------------//

/*!
 * \brief Integrates the tensor (Omega Omega) over the unit sphere. 
 *
 * The solution to this integral is a tensor vector with ndims^2 elements
 * The off-diagonal elements should be zero.  The diagonal elements
 * should have the value sumwt/3.  
 *
 * We actually return a 1D vector whose length is ndims^2.  The "diagonal"
 * elements are 0, 4 and 8.
 *
 * The integral is actually calculated as a quadrature sum over all
 * directions of the quantity:
 *
 *     Omega(m) Omega(m) wt(m)
 *
 * The quantity ( Omega Omega ) is a tensor quantity.
 */
vector<double> Quadrature::iOmegaOmegaDomega() const {
    size_t ndims = dimensionality();
    // The size of the solution tensor will be ndims^2.
    // The solution is returned as a vector and not a tensor.  The diagonal
    // elements of the tensor are elements 0, 4 and 8 of the vector.
    vector<double> integral( ndims*ndims );
    // initialize the solution to zero.
    for ( size_t i = 0; i < ndims*ndims; ++i) 
	integral[0] = 0.0;
    // We are careful to only compute the terms of the tensor solution that
    // are available for the current dimensionality of the quadrature set.
    switch (ndims) {
    case 1:
	for ( size_t i = 0; i < getNumAngles(); ++i ) {
	    integral[0] += wt[i]*mu[i]*mu[i];
	}
	break;
    case 2:
	for ( size_t i = 0; i < getNumAngles(); ++i ) {
	    integral[0] += wt[i]*mu[i]*mu[i];
	    integral[1] += wt[i]*mu[i]*eta[i];
	    integral[2] += wt[i]*eta[i]*mu[i];
	    integral[3] += wt[i]*eta[i]*eta[i];
	}
	break;
    case 3:
	for ( size_t i = 0; i < getNumAngles(); ++i ) {
	    integral[0] += wt[i]*mu[i]*mu[i];
	    integral[1] += wt[i]*mu[i]*eta[i];
	    integral[2] += wt[i]*mu[i]*xi[i];
	    integral[3] += wt[i]*eta[i]*mu[i];
	    integral[4] += wt[i]*eta[i]*eta[i];
	    integral[5] += wt[i]*eta[i]*xi[i];
	    integral[6] += wt[i]*xi[i]*mu[i];
	    integral[7] += wt[i]*xi[i]*eta[i];
	    integral[8] += wt[i]*xi[i]*xi[i];
	}
	break;
    default:
	Insist(false,"Number of spatial dimensions must be 1, 2 or 3.");
    }
    return integral;
}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of Quadrature.cc
//---------------------------------------------------------------------------//
