//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q3DTriChebyshevLegendre.cc
 * \author James Warsa
 * \date   Fri Jun  9 13:52:25 2006
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "Q1DGaussLeg.hh"
#include "Q3DTriChebyshevLegendre.hh"
#include "Ordinate.hh"

namespace rtt_quadrature
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructs a 3D Tri Chebyshev Legendre quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder*snOrder)
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2*PI).
 */
//---------------------------------------------------------------------------//
Q3DTriChebyshevLegendre::Q3DTriChebyshevLegendre( size_t sn_order_, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( sn_order_, norm_, qm_ ),
      numOrdinates ( sn_order_ * (sn_order_+2) )
{
    using std::fabs;
    using std::sqrt;
    using std::cos;
    using rtt_dsxx::soft_equiv;
    
    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    Insist( snOrder%2 == 0, "SN order must be even");

    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates, 0.0);
    eta.resize(numOrdinates, 0.0);
    xi.resize(numOrdinates, 0.0);
    wt.resize(numOrdinates, 0.0);

    Q1DGaussLeg gauss(snOrder, 2.0, interpModel);    

    // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)

    unsigned icount=0;

    for (unsigned i=0; i<snOrder/2; ++i)
    {
        double xmu=gauss.getMu(i);
        double xwt=gauss.getWt(i);
        double xsr=sqrt(1.0-xmu*xmu);

        unsigned const k=2*(i+1);
            
        for (unsigned j=0; j<k/2; ++j)
        {
            unsigned ordinate=icount;
                
            mu[ordinate] = xmu;
            eta[ordinate] = xsr*cos(rtt_units::PI*(2.0*j+1.0)/k/2.0);
            wt[ordinate] = xwt/k;

            ++icount;
        }
    }

    // The number of quadrature levels is equal to the requested SN order.
    size_t octantOrdinates( numOrdinates/8 );  // 8 octants in 3D.

    // Should only have computed one octant
    Insist(icount == octantOrdinates, "Computed an unexpected number of ordinates for the first octant");

    // Evaluate mu and eta for octants 2-4
    for(size_t octant=2; octant<=4; ++octant)
    {
	for(size_t n=0; n<=octantOrdinates-1; ++n) 
	{
            ++icount;
	    unsigned const m = (octant-1)*octantOrdinates+n;
	    switch (octant) {
	    case 2:
		mu[m]  = -mu[n];
		eta[m] =  eta[n];
		wt[m]  =  wt[n];
		break;
		
	    case 3:
		mu[m]  = -mu[n];
		eta[m] = -eta[n];
		wt[m]  =  wt[n];
		break;
		
	    case 4:
		mu[m]  =  mu[n];
		eta[m] = -eta[n];
		wt[m]  =  wt[n];
		break;
	    default:
		Insist(false,"Octant value should only be 2, 3 or 4 in this loop.");
		break;
	    }
	}
    }
    
    // Evaluate mu and eta for octants 5-8
    for( size_t n=0; n<=numOrdinates/2-1; ++n)
    {
        ++icount;
	mu[n+numOrdinates/2]  = mu[n];
	eta[n+numOrdinates/2] = eta[n];
	wt[n+numOrdinates/2]  = wt[n];
    }
    
    // Evaluate xi for all octants
    for( size_t n=0;n<=numOrdinates/2-1;++n)
	xi[n] = sqrt(1.0-mu[n]*mu[n]-eta[n]*eta[n]);
    
    for(size_t n=0;n<=numOrdinates/2-1;++n)
	xi[n+numOrdinates/2] = -xi[n];

    Insist(icount == numOrdinates, "Computed an unexpected number of ordinates");

    // Normalize the quadrature set
    double wsum = 0.0;
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wsum = wsum + wt[ordinate];
    
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wt[ordinate] = wt[ordinate]*(norm/wsum);

    // Verify that the quadrature meets our integration requirements.
    Ensure( soft_equiv(iDomega(),norm) );

    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    Ensure( soft_equiv(iod[1],0.0) );
    Ensure( soft_equiv(iod[2],0.0) );
		    
    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( soft_equiv(iood[0],norm/3.0) );  // mu*mu
    Ensure( soft_equiv(iood[1],0.0) ); // mu*eta
    Ensure( soft_equiv(iood[2],0.0) ); // mu*xi
    Ensure( soft_equiv(iood[3],0.0) ); // eta*mu
    Ensure( soft_equiv(iood[4],norm/3.0) ); // eta*eta
    Ensure( soft_equiv(iood[5],0.0) ); // eta*xi
    Ensure( soft_equiv(iood[6],0.0) ); // xi*mu
    Ensure( soft_equiv(iood[7],0.0) ); // xi*eta
    Ensure( soft_equiv(iood[8],norm/3.0) ); // xi*xi

    // Copy quadrature data { mu, eta } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims = dimensionality();
    for ( size_t ordinate = 0; ordinate < numOrdinates; ++ordinate )
    {
	omega[ordinate].resize(ndims);
	omega[ordinate][0] = mu[ordinate];
	omega[ordinate][1] = eta[ordinate];
	omega[ordinate][2] = xi[ordinate];
    }

    //display();

} // end of Q3DTriChebyshevLegendre constructor

//---------------------------------------------------------------------------//

void Q3DTriChebyshevLegendre::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;
    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t    eta       \t    xi        \t     wt      " << endl;
    cout << "  --- \t------------- \t------------- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ix = 0; ix < mu.size(); ++ix ) {
	cout << "   "
	     << ix << "\t"
	     << setprecision(10) << mu[ix]  << "\t"
	     << setprecision(10) << eta[ix] << "\t"
	     << setprecision(10) << xi[ix]  << "\t"
	     << setprecision(10) << wt[ix]  << endl;
	sum_wt += wt[ix];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q3DTriChebyshevLegendre.cc
//---------------------------------------------------------------------------//
