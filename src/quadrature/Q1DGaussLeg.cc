//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1DGaussLeg.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 09:35:03 2004
 * \brief  1D Gauss Legendre Quadrature
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>

#include "units/PhysicalConstants.hh"
#include "Quadrature.hh"
#include "Q1DGaussLeg.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 1D Gauss Legendre quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed. 
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2.0).
 */

Q1DGaussLeg::Q1DGaussLeg( size_t n, double norm_ ) 
    : Quadrature( n, norm_ ), numAngles( n )
{
    // We require the sn_order to be greater than zero.
    Require( n > 0 );
    // We require the normalization constant to be greater than zero.
    Require( norm > 0.0 );

    // size the member data vectors
    mu.resize(n);
    wt.resize(n);

    double mu1 = -1.0;
    double mu2 =  1.0;

    size_t m = (n+1)/2;  // number of angles in half range.
    double mu_m = 0.5 * (mu2+mu1);
    double mu_l = 0.5 * (mu2-mu1);
    
    for ( size_t i = 1; i <= m; ++i ) {
	double z1, pp, z = cos( rtt_units::PI * ( i-0.25 ) / ( n+0.5 ));
	do {
	    double p1 = 1.0, p2 = 0.0;
	    for ( size_t j=1; j<=n; j++ ) {
		double p3 = p2;
		p2 = p1;
		p1 = (( 2.0*j - 1.0 ) * z * p2 - ( j - 1.0 ) * p3 ) / j;
	    }
	    pp = n * ( z*p1 - p2 ) / ( z * z - 1.0 );
	    z1 = z;
	    z = z1-p1/pp;
	} while ( fabs( z-z1 ) > TOL );

	// fill the vectors
	mu[ i-1 ] = mu_m - mu_l * z;
	wt[ i-1 ] = 2.0 * mu_l / (( 1.0 - z*z ) * pp*pp );
	mu[ n-i ] = mu_m + mu_l * z;
	wt[ n-i ] = wt[ i-1 ];
    }	

    // Sanity Checks: 
    // Verify that the 0th, 1st and 2nd moments integrate to 2, (0,0,0) and
    // 2/3*((1,0,0),(0,1,0),(0,0,1)) respectively.

    double sumwt = 0.0;
    //    double gamma = 0.0;
    //    double zeta  = 0.0;

    for ( size_t i = 0; i < n; ++i )
    {
	sumwt += wt[i];
	//	gamma += wt[i] * mu[i];
	//	zeta  += wt[i] * mu[i] * mu[i];
    }

    // The quadrature weights should sum to 2.0
    Ensure( fabs(iDomega()-2.0) <= TOL );
    // The integral of mu over all angles should be zero.
    Ensure( fabs(iOmegaDomega()[0]) <= TOL );
    // The integral of mu^2 should be 2/3.
    Ensure( fabs(iOmegaOmegaDomega()[0]-2.0/3.0) <= TOL );

    // If norm != 2.0 then renormalize the weights to the required values. 
    if ( fabs(norm-2.0) >= TOL ) 
    {
	double c = norm/sumwt;
	for ( size_t i=0; i < numAngles; ++i )
	    wt[i] = c * wt[i];
    }
    
    // make a copy of the data into the omega vector < vector < double > >
    omega.resize( n );
    for ( size_t i=0; i<n; ++i )
    {
	// This is a 1D set.
	omega[i].resize(1);
	omega[i][0] = mu[i];
    }

} // end of Q1DGaussLeg() constructor.

//---------------------------------------------------------------------------//

void Q1DGaussLeg::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;	

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t     wt      " << endl;
    cout << "  --- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ix = 0; ix < mu.size(); ++ix ) {
	cout << "   "
	     << setprecision(5)  << ix     << "\t"
	     << setprecision(10) << mu[ix] << "\t"
	     << setprecision(10) << wt[ix] << endl;
	sum_wt += wt[ix];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}


} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q1DGaussLeg.cc
//---------------------------------------------------------------------------//
