//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1DDoubleGauss.cc
 * \author James Warsa
 * \date   Fri Sep 16 15:45:26 2005
 * \brief  1D Double-Gauss Quadrature
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "Quadrature.hh"
#include "gauleg.hh"
#include "Q1DDoubleGauss.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 1D double-Gauss quadrature object.
 *
 * Calculates a double-Gauss quadrature using abscissas and weights 
 *
 * For normalized lower and upper limits of integration -1.0 & 1.0, and given
 * n, this routine first calculates arrays of length n/2 containing abscissas and
 * weights of a Gauss-Legendre (n/2)-point quadrature formula and maps that
 * onto the two half-ranges [-1,0] and [0, 1]. 
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed. 
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2.0).
 */

Q1DDoubleGauss::Q1DDoubleGauss( size_t numGaussPoints, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( numGaussPoints, norm_, qm_ ), 
      numOrdinates( numGaussPoints )
{
    using rtt_dsxx::soft_equiv;

    // We require the sn_order to be greater than zero.
    Require( numGaussPoints > 0 );
    // And that it is even
    Require( numGaussPoints%2 == 0 );
    // We require the normalization constant to be greater than zero.
    Require( norm > 0.0 );

    unsigned const n(numGaussPoints);
    unsigned const n2(n/2);

    // size the data vectors

    mu.resize(numGaussPoints);
    wt.resize(numGaussPoints);

    if (n2 == 1)  
    {
        // 2-point double Gauss is just Gauss

        double const mu1(-1); // range of direction
        double const mu2(1);
        gauleg( mu1, mu2, mu, wt, numGaussPoints );
        
    }
    else
    {
        // Create an N/2-point Gauss quadrature on [-1,1] 

        Check( n2%2 == 0 );

        double const mu1(-1); // range of direction
        double const mu2(1);
        std::vector< double > muH;
        std::vector< double > wtH;
        gauleg( mu1, mu2, muH, wtH, n2 );
        
        // map the quadrature onto the two half-ranges
        
        for (unsigned m=0; m<n2; ++m)
        {
            // Map onto [-1,0] then skew-symmetrize
            // (ensuring ascending order on [-1, 1])
            
            mu[m] = 0.5*(muH[m] - 1.0);
            wt[m] = 0.5*wtH[m];
            
            mu[n-m-1] = -mu[m];
            wt[n-m-1] =  wt[m];
        }
    }

    // Compute and store the sum of the weights
    
    double sumwt( std::accumulate( wt.begin(), wt.end(), // range
                                   0.0 ) );              // init value.

    // Sanity Checks: always none at present

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

} // end of Q1DDoubleGauss() constructor.

//---------------------------------------------------------------------------//

void Q1DDoubleGauss::display() const 
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
//                 end of Q1DDoubleGauss.cc
//---------------------------------------------------------------------------//
