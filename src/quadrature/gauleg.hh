//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/gauleg.hh
 * \author Kent Budge
 * \date   Tue Sep 14 13:16:09 2004
 * \brief  Gauss-Legendre quadrature
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_quadrature_gauleg_hh
#define rtt_quadrature_gauleg_hh

#include <limits>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_quadrature
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Gauss-Legendre quadrature
 *
 * Calculate abscissae and weights for Gauss-Legendre quadrature:
 * \f[
 * \int_{x_1}^{x_2}{f(x) dx} = \sum\limits_{j=0}^{N - 1}w_j f(x_j)
 * \f]
 *
 * We will use Newton's method of root finding to determine the
 * abscissas. This algorithm is a modified form of that found in "Numerical
 * Recipes in C."
 *
 * The abcissas are the roots of this recurrence relation:
 * \f[
 * (j+1)P_{j+1} = (2j-1)xP_j - jP_{j-1}
 * \f]
 *
 * The weights are determined from the relation:
 * \f[
 * w_j = \frac{2}{(1-x^2_j)[P'_N(x_j)]^2}
 * \f]
 *
 * This routine scales the range of integration from \f$(x_1,x_2)\f$ to (-1,1)
 * and provides the abscissas \f$ x_j\f$ and weights \f$ w_j \f$ for the
 * Gaussian formula provided above.
 * 
 * \arg \a FieldVector A random access container on a field type.
 *
 * \param x1 Start of integration interval
 * \param x2 End of integration interval
 * \param x On return, contains abscissae \f$x_j\f$ for quadrature.
 * \param w On return, contains weights \f$w_j\f$ for quadrature.
 * \param N Number of points in quadrature.
 * 
 */
template< typename FieldVector >
void gauleg( double const x1, // expect FieldVector::value_type to be promoted to double.
             double const x2, // expect FieldVector::value_type to be promoted to double.
             FieldVector &x,
             FieldVector &w,
             unsigned const n )
{
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;
    using std::cos;
    using std::numeric_limits;

    typedef typename FieldVector::value_type Field;

    Require( n > 0 );
    Require( x2 > x1 );

    // convergence tolerance
    Field const tolerance( 100*numeric_limits< Field >::epsilon() );
    
    x.resize(n);
    w.resize(n);

    // number of Gauss points in the half range.
    // The roots are symmetric in the interval.  We only need to search for
    // half of them.
    unsigned const numHrGaussPoints( (n+1)/2 );

    // mid-point of integration range
    Field const mu_m( 0.5*(x2+x1) );
    // length of half the integration range.
    Field const mu_l( 0.5*(x2-x1) );
    
    // Loop over the desired roots.
    for ( size_t iroot=0; iroot<numHrGaussPoints; ++iroot)
    {
	// Approximate the i-th root.
	Field z( cos( PI * ( iroot+0.75 ) / ( n+0.5 )) );

        // Temporary storage;
        Field z1, pp;
        
	do // Use Newton's method to refine the value for the i-th root.  
	{
	    Field p1( 1 );
	    Field p2( 0 );

	    // This loop represents the recurrence relation needed to obtain
	    // the Legendre polynomial evaluated at z.
	    for( unsigned j=0; j<n; ++j )
	    {
		Field const p3 = p2;
		p2 = p1;
		p1 = (( 2*j+1) * z * p2 - j*p3 ) / (j+1);
	    }
	    
	    // p1 is now the desired Legendre polynomial evaluated at z. We
	    // next compute pp, its derivative, by a standard relation
	    // involving also p2, the polynomial of one lower order.
	    pp = n * (z*p1-p2)/(z*z-1.0);
	    z1 = z;
	    
	    // Newton's Method
	    z = z1 - p1/pp;

	} while( ! soft_equiv(z,z1, tolerance ) );

	// Roots will be between -1 and 1.0 and symmetric about the origin. 
        int const idxSymPart( n - iroot - 1 );

        // Now, scale the root to tthe desired interval and put in its
        // symmetric counterpart.
	x[ iroot ]      = mu_m - mu_l*z;       
	x[ idxSymPart ] = mu_m + mu_l*z;       

	// Compute the associated weight and its symmetric counterpart.
        w[ iroot ]      = 2*mu_l / ((1-z*z)*pp*pp); 
        w[ idxSymPart ] = w[ iroot ];
    }	
    
// Kent's original code:
    
//     unsigned const m = (n+1)/2;
//     double const xm = 0.5*(x2+x1);
//     double const xl = 0.5*(x2-x1);
//     for (unsigned i=0; i<m; ++i)
//     {
//         double z = cos(PI*(i+0.75)/(n+0.5));
//         double z1, pp;
//         do
//         {
//             double p1 = 1.0;
//             double p2 = 0.0;
//             for (unsigned j=0; j<n; ++j)
//             {
//                 double const p3 = p2;
//                 p2 = p1;
//                 p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1);
//             }
//             pp = n*(z*p1-p2)/(z*z-1.0);
//             z1 = z;
//             z = z1-p1/pp;
//         }
//         while (fabs(z-z1)>100*numeric_limits<double>::epsilon());
//         x[i] = xm-xl*z;
//         x[n-1-i] = xm+xl*z;
//         w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
//         w[n-1-i] = w[i];
//     }

    return;
}

} // end namespace rtt_quadrature

#endif // rtt_quadrature_gauleg_hh

//---------------------------------------------------------------------------//
//              end of quadrature/gauleg.hh
//---------------------------------------------------------------------------//
