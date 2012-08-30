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

#include <cmath>
#include <sstream>

#include "Quadrature.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Integrates dOmega over the unit sphere. (The sum of quadrature weights.)
 */
double Quadrature::iDomega() const {
    double integral = 0.0;
    size_t N( wt.size() );
    for ( size_t i=0; i<N; ++i )
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
    vector<double> integral(3,0.0);
    Require( wt.size() > 0 );
    Require( mu.size() == wt.size() );
    unsigned n(mu.size());
    for( unsigned i=0; i<n; ++i )
        integral[0] += wt[i]*mu[i];
    n = eta.size();
    Check( n <= wt.size() );
    for( unsigned i=0; i<n; ++i )
        integral[1] += wt[i]*eta[i];
    n = xi.size();
    Check( n <= wt.size() );
    for( unsigned i=0; i<n; ++i )
        integral[2] += wt[i]*xi[i];    
//     switch( ndims ) {
//     case 3:
// 	for ( size_t i = 0; i < getNumOrdinates(); ++i )
// 	    integral[2] += wt[i]*xi[i];
// 	//lint -fallthrough
//     case 2:
// 	for ( size_t i = 0; i < getNumOrdinates(); ++i )
// 	    integral[1] += wt[i]*eta[i];
// 	//lint -fallthrough
//     case 1:
// 	for ( size_t i = 0; i < getNumOrdinates(); ++i )
// 	    integral[0] += wt[i]*mu[i];
// 	break;
//     default:
// 	Insist(false,"Number of spatial dimensions must be 1, 2 or 3.");
//    }
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
//     // may want full tensor for 2D quadrature set.
//     if( eta.size() > 0 && eta.size() == xi.size() )
//         ndims = 3;
    // The size of the solution tensor will be ndims^2.
    // The solution is returned as a vector and not a tensor.  The diagonal
    // elements of the tensor are elements 0, 4 and 8 of the vector.
    vector<double> integral( ndims*ndims, 0.0 );

    // We are careful to only compute the terms of the tensor solution that
    // are available for the current dimensionality of the quadrature set.
    if( xi.size()>0 && eta.size()>0 )
    {
        for ( size_t i = 0; i < getNumOrdinates(); ++i )
        {
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
    }
    else if( xi.size()>0 )
    {
        for ( size_t i = 0; i < getNumOrdinates(); ++i )
        {
            integral[0] += wt[i]*mu[i]*mu[i];
            integral[1] += wt[i]*mu[i]*xi[i];
            integral[2] += wt[i]*xi[i]*mu[i];
            integral[3] += wt[i]*xi[i]*xi[i];
        }
    }
    else if( eta.size()>0 )
    {
        for ( size_t i = 0; i < getNumOrdinates(); ++i )
        {
            integral[0] += wt[i]*mu[i]*mu[i];
            integral[1] += wt[i]*mu[i]*eta[i];
            integral[2] += wt[i]*eta[i]*mu[i];
            integral[3] += wt[i]*eta[i]*eta[i];
        }
    }
    else
    {
	for ( size_t i = 0; i < getNumOrdinates(); ++i )
	    integral[0] += wt[i]*mu[i]*mu[i];
    }
    return integral;
}

//---------------------------------------------------------------------------//
void Quadrature::renormalize(const double new_norm)
{
    Require(new_norm > 0);
    Require(norm > 0);

    double const c=new_norm/norm;
    
    for ( size_t i = 0; i < getNumOrdinates(); ++i )
    {
        wt[i] = c * wt[i];
    }

    // re-set normalization value
    norm = new_norm;
}

//---------------------------------------------------------------------------//
string Quadrature::as_text(string const &indent) const
{
    using namespace std;

    stringstream text;
    text << indent + "  type, " << parse_name()
         << indent + "  order, " << snOrder
         << indent + "end";

    return text.str();
}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                         end of Quadrature.cc
//---------------------------------------------------------------------------//
