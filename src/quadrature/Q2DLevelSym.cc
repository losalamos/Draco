//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q2DLevelSym.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:27:00 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>

#include "Q3DLevelSym.hh"
#include "Q2DLevelSym.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructs a 2D Level Symmetric quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of angles = (snOrder+2)*snOrder/2.
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2*PI).
 */
//---------------------------------------------------------------------------//
Q2DLevelSym::Q2DLevelSym( size_t sn_order_, double norm_ ) 
    : Quadrature( sn_order_, norm_ ), numAngles (sn_order_ * (sn_order_+2)/2)
{
    using std::fabs;

    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    // Insist ( snOrder%2 == 0, "LS Quad must have an even SN order." );
    Require ( snOrder%2 == 0 );
    // Insist ( snOrder >= 2 && snOrder <= 24, "LS Quad must have a SN order between 2 and 24." );
    Require ( snOrder >= 2 && snOrder <= 24 );

    // Force the direction vectors to be the correct length.
    mu.resize(numAngles);
    eta.resize(numAngles);
    wt.resize(numAngles);

    // Use the 3D level-symmetric quadrature to construct the 2D one.
    // In the future we probably should reverse the dependence - we should
    // use the 2D quadrature to construct the 3D one, rather than depending
    // on a particular data layout in the 3D quadrature.
    Q3DLevelSym quad3D(snOrder, 2.0*norm);
    for( size_t angle = 0; angle < numAngles; ++angle )
    {
	mu[angle] = quad3D.getMu(angle);
	eta[angle] = quad3D.getEta(angle);
	wt[angle] = quad3D.getWt(angle);
    }

    // Normalize the quadrature set
    double wsum = 0.0;
    for(size_t angle = 0; angle < numAngles; ++angle)
	wsum = wsum + wt[angle];
    
    for(size_t angle = 0; angle < numAngles; ++angle)
	wt[angle] = wt[angle]*(norm/wsum);

    // Verify that the quadrature meets our integration requirements.
    Ensure( fabs(iDomega()-norm) <= TOL );

    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( fabs(iod[0]) <= TOL );
    Ensure( fabs(iod[1]) <= TOL );

    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( fabs(iood[0]-norm/3.0 ) <= TOL );  // mu*mu
    Ensure( fabs(iood[1]) <= TOL ); // mu*eta
    Ensure( fabs(iood[2]) <= TOL ); // eta*mu
    Ensure( fabs(iood[3]-norm/3.0) <= TOL ); // eta*eta

    // Copy quadrature data { mu, eta } into the vector omega.
    omega.resize( numAngles );
    size_t ndims = dimensionality();
    for ( size_t angle = 0; angle < numAngles; ++angle )
    {
	omega[angle].resize(ndims);
	omega[angle][0] = mu[angle];
	omega[angle][1] = eta[angle];
    }
} // end of Q2DLevelSym() constructor.

//---------------------------------------------------------------------------//

void Q2DLevelSym::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t    eta       \t     wt      " << endl;
    cout << "  --- \t------------- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t angle = 0; angle < mu.size(); ++angle )
    {
	cout << "   "
	     << angle << "\t"
	     << setprecision(10) << mu[angle]  << "\t"
	     << setprecision(10) << eta[angle] << "\t"
	     << setprecision(10) << wt[angle]  << endl;
	sum_wt += wt[angle];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q2DLevelSym.cc
//---------------------------------------------------------------------------//
