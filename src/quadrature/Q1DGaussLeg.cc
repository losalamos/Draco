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
#include <cmath>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "Quadrature.hh"
#include "Q1DGaussLeg.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 1D Gauss Legendre quadrature object.
 *
 * Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian
 * Quadrature integration of polynomial functions.
 *
 * For normalized lower and upper limits of integration -1.0 & 1.0, and given
 * n, this routine calculates, arrays xabsc(1:n) and weig(1:n) of length n,
 * containing the abscissas and weights of the Gauss-Legendre n-point
 * quadrature formula.  For detailed explanations finding weights &
 * abscissas, see "Numerical Recipes in Fortran.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed. 
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2.0).
 */

Q1DGaussLeg::Q1DGaussLeg( size_t numGaussPoints, double norm_ ) 
    : Quadrature( numGaussPoints, norm_ ), 
      numAngles( numGaussPoints )
{
    using rtt_dsxx::soft_equiv;

    // We require the sn_order to be greater than zero.
    Require( numGaussPoints > 0 );
    // We require the normalization constant to be greater than zero.
    Require( norm > 0.0 );

    // size the member data vectors
    mu.resize(numGaussPoints);
    wt.resize(numGaussPoints);

    double mu1 = -1.0;  // minimum value for mu
    double mu2 =  1.0;  // maximum value for mu

    // number of Gauss points in the half range.
    // The roots are symmetric in the interval.  We only need to search for
    // half of them.
    unsigned const numHrGaussPoints( (numGaussPoints+1)/2 );

    double mu_m = 0.5 * (mu2+mu1);
    double mu_l = 0.5 * (mu2-mu1);
    double z1,pp;

    // Loop over the desired roots.
    for ( size_t iroot=0; iroot<numHrGaussPoints; )
    {
	++iroot;

	// Approximate the i-th root.
	double z( std::cos( rtt_units::PI * ( iroot-0.25 ) 
		       / ( numGaussPoints+0.5 )) );

	do // Use Newton's method to refine the value for the i-th root.  
	{
	    double p1( 1.0 );
	    double p2( 0.0 );
	    double p3;

	    // This loop represents the recurrence relation needed to obtain
	    // the Legendre polynomial evaluated at z.
	    for( unsigned j=0; j<numGaussPoints; )
	    {
		++j;
		p3 = p2;
		p2 = p1;
		p1 = (( 2.0*j-1.0) * z * p2 - (j-1) * p3 ) / j;
	    }
	    
	    // p1 is now the desired Legendre polynomial evaluated at z. We
	    // next compute pp, its derivative, by a standard relation
	    // involving also p2, the polynomial of one lower order.
	    pp = numGaussPoints * (z*p1-p2)/(z*z-1.0);
	    z1 = z;
	    
	    // Newton's Method
	    z = z1 - p1/pp;   

	} while( ! soft_equiv(z,z1) );

	// Roots wil be between -1 and 1.0 and symmetric about the origin. 
	mu[ iroot-1 ]              = -z;       
	mu[ numGaussPoints-iroot ] = z;       

	// Compute the associated weight and its symmetric counterpart.
        wt[ iroot-1 ]              = 2.0/((1.0-z*z)*pp*pp); 
        wt[ numGaussPoints-iroot ] = wt[iroot-1];
    
    }	

    // Sanity Checks: 

    // Verify that the 0th, 1st and 2nd moments integrate to 2, (0,0,0) and 
    // 2/3*((1,0,0),(0,1,0),(0,0,1)) respectively.

    double sumwt = 0.0;
    //    double gamma = 0.0;
    //    double zeta  = 0.0;

    for ( size_t i = 0; i < numGaussPoints; ++i )
    {
	sumwt += wt[i];
	//	gamma += wt[i] * mu[i];
	//	zeta  += wt[i] * mu[i] * mu[i];
    }

    // The quadrature weights should sum to 2.0
    Ensure( soft_equiv(iDomega(),2.0) );
    // The integral of mu over all angles should be zero.
    Ensure( soft_equiv(iOmegaDomega()[0],0.0) );
    // The integral of mu^2 should be 2/3.
    Ensure( soft_equiv(iOmegaOmegaDomega()[0],2.0/3.0) );

    // If norm != 2.0 then renormalize the weights to the required values. 
    if( soft_equiv(norm,2.0) ) 
    {
	double c = norm/sumwt;
	for ( size_t i=0; i < numAngles; ++i )
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
