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
#include <numeric>

#include "ds++/Soft_Equivalence.hh"

#include "Quadrature.hh"
#include "Q1DGaussLeg.hh"
#include "gauleg.hh"

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

Q1DGaussLeg::Q1DGaussLeg( size_t numGaussPoints, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( numGaussPoints, norm_, qm_ ), 
      numOrdinates( numGaussPoints )
{
    using rtt_dsxx::soft_equiv;

    // We require the sn_order to be greater than zero.
    Require( numGaussPoints > 0 );
    // We require the normalization constant to be greater than zero.
    Require( norm > 0.0 );

    double const mu1(-1); // range of direction
    double const mu2( 1);
    gauleg( mu1, mu2, mu, wt, numGaussPoints );
    
    double sumwt( std::accumulate( wt.begin(), wt.end(), // range
                                   0.0 ) );              // init value.

    // The quadrature weights should sum to 2.0
    Ensure( soft_equiv(iDomega(),2.0) );
    // The integral of mu over all ordinates should be zero.
    Ensure( soft_equiv(iOmegaDomega()[0],0.0) );
    // The integral of mu^2 should be 2/3.
    Ensure( soft_equiv(iOmegaOmegaDomega()[0],2.0/3.0) );

    // If norm != 2.0 then renormalize the weights to the required values. 
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

    return;
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
