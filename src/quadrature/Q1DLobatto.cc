//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /Q1DLobatto.cc
 * \author James Warsa
 * \date   Fri Sep  2 10:30:02 2005
 * \brief  Creates one-dimensional Lobatto quadrature.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <limits>

#include <gsl/gsl_sf_legendre.h>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "Quadrature.hh"
#include "Q1DLobatto.hh"

namespace rtt_quadrature
{
/*!
 * \brief Constructs a 1D Lobatto quadrature object.
 *
 * Calculation of Lobatto abscissas and weights for 
 * Quadrature integration of polynomial functions.
 *
 * For normalized lower and upper limits of integration -1.0 & 1.0, and given
 * n, this routine calculates, arrays xabsc(1:n) and weig(1:n) of length n,
 * containing the abscissas and weights of the Lobatto n-point
 * quadrature formula.  For detailed explanations finding weights &
 * abscissas, see "Numerical Recipes in Fortran.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed. 
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2.0).
 */

Q1DLobatto::Q1DLobatto( size_t numGaussPoints, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( numGaussPoints, norm_, qm_ ), 
      numOrdinates( numGaussPoints )
{
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;
    using std::cos;
    using std::numeric_limits;

    // We require the sn_order to be greater than zero.
    Require( numGaussPoints > 0 );
    // And that it is even
    Require( numGaussPoints%2 == 0 );
    // We require the normalization constant to be greater than zero.
    Require( norm > 0.0 );

    double const tolerance( 100*std::numeric_limits< double >::epsilon() );

    // size the member data vectors
    unsigned const N( numGaussPoints );
    mu.resize( N );
    wt.resize( N );

    double const mu1( -1.0 );  // minimum value for mu
    double const mu2(  1.0 );  // maximum value for mu

    // Number of Gauss points in the half range.
    // The roots are symmetric in the interval.  We only need to search for
    // half of them.
    unsigned const numHrGaussPoints( numGaussPoints/2 );

    mu[0]   = mu1;
    mu[N-1] = mu2;

    // Loop over the desired roots.
    for ( size_t iroot=0; iroot<numHrGaussPoints-1; ++iroot)
    {
	// Approximate the i-th root.
	double z( cos( PI * ( iroot-0.25 ) / ( (N-2)+0.5 )) );
        double z1;

	do // Use Newton's method to refine the value for the i-th root.  
	{
            // P_{N-1}(z)
            double const pnm1 ( gsl_sf_legendre_Pl(N-1, z) );
            // P_{N}(z)
            double const pn   ( gsl_sf_legendre_Pl(N  , z) );
            // P_{N+1}(z)
            double const pnp1 ( gsl_sf_legendre_Pl(N+1, z) );

            // dP/dz _{N-1}(z)
            double const   pp( (N)   * (z*pnm1 - pn) / (1.0 - z*z) );
            // dP/dz _{N}(z)
            double const  pp1( (N+1) * (z*pn - pnp1) / (1.0 - z*z) );

            // d2P/dz2 _{N}(z)
            double const  pdp( ( (N) * (z*pp + pnm1 - pp1) + 2*z*pp ) / (1.0 - z*z) );

	    // update 

	    z1 = z;
	    z = z1 - pp/pdp;   

 	} while( ! soft_equiv(z,z1, tolerance) );

	// Roots will be in [-1,1], symmetric about the origin. 
	mu[ iroot+1 ]                 = -z;       
	mu[ numGaussPoints-iroot-2] =  z;       
    }	

    // Loop over the quadrature points to compute weights.
    for ( size_t m=0; m<numHrGaussPoints; ++m)
    {
        double const z(mu[m]);
        double const p( gsl_sf_legendre_Pl(N-1, z) );

	// Compute the associated weight and its symmetric counterpart.
        wt[ m ]                 = 2.0/N/(N-1)/p/p;
        wt[ numGaussPoints-m-1] = wt[m];
    }

    // Sanity Checks: none at present

    double sumwt = 0.0;
    for ( size_t i = 0; i < numGaussPoints; ++i )
	sumwt += wt[i];

    if( !soft_equiv(norm,2.0) ) 
    {
	double c = norm/sumwt;
	for ( size_t i=0; i < numOrdinates; ++i )
	    wt[i] = c * wt[i];
    }
    
    // make a copy of the data into the omega vector < vector < double > >
    omega.resize( numGaussPoints );
    for ( size_t i=0; i<numGaussPoints; ++i )
    {
	// This is a 1D set.
	omega[i].resize(1);
	omega[i][0] = mu[i];
    }

} // end of Q1DLobatto() constructor.

//---------------------------------------------------------------------------//

void Q1DLobatto::display() const 
{
    using namespace std;

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "     m           mu             wt     " << endl;
    cout << "  -------  -------------  -------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ix = 0; ix < mu.size(); ++ix ) {
        cout.setf(ios::right);
        cout.setf(ios::fixed); 
        cout.setf(ios::showpoint);
	cout << setw(7)  << setprecision(5) 
             << "  " << ix;
	cout << setprecision(10) 
             << "   " << mu[ix] 
             << "   " << wt[ix] 
             << endl;
        cout.unsetf(ios::right);
        cout.unsetf(ios::showpoint);
        cout.unsetf(ios::floatfield);
	sum_wt += wt[ix];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q1DLobatto.cc
//---------------------------------------------------------------------------//
