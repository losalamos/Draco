//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q2DTriChebyshevLegendre.cc
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
#include "Q2DTriChebyshevLegendre.hh"
#include "Ordinate.hh"

namespace rtt_quadrature
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructs a 2D Tri Chebyshev Legendre quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder*snOrder)
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2*PI).
 */
//---------------------------------------------------------------------------//
Q2DTriChebyshevLegendre::Q2DTriChebyshevLegendre( size_t sn_order_, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( sn_order_, norm_, qm_ ), numOrdinates (sn_order_ * (sn_order_ + 2) / 2)
{
    using std::fabs;
    using std::sqrt;
    using std::cos;
    using rtt_dsxx::soft_equiv;
    
    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    Insist( snOrder%2 == 0, "SN order must be even");

    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates);
    xi.resize(numOrdinates);
    wt.resize(numOrdinates);

    Q1DGaussLeg gauss(snOrder, 2.0, interpModel);    

    // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)

    unsigned icount=0;

    for (unsigned i=0; i<snOrder/2; ++i)
    {
        double xmu=gauss.getMu(i);
        double xwt=gauss.getWt(i);
        double xsr=sqrt(1.0-xmu*xmu);

        unsigned const k=2*(i+1);
            
        for (unsigned j=0; j<k; ++j)
        {
            unsigned ordinate=icount;
                
            xi[ordinate] = xmu;
            mu[ordinate]  = xsr*cos(rtt_units::PI*(2.0*j+1.0)/k/2.0);
            wt[ordinate]  = xwt/k;

            ++icount;
        }

        unsigned const ii=snOrder-i-1;

        xmu=gauss.getMu(ii);
        xwt=gauss.getWt(ii);
        xsr=sqrt(1.0-xmu*xmu);

        for (unsigned j=0; j<k; ++j)
        {
            unsigned ordinate=icount;
                
            xi[ordinate] = xmu;
            mu[ordinate]  = xsr*cos(rtt_units::PI*(2.0*j+1.0)/k/2.0);
            wt[ordinate]  = xwt/k;

            ++icount;
        }
    }

    Insist(icount == numOrdinates, "Computed an unexpected number of ordinates");

    // Normalize the quadrature set
    double wsum = 0.0;
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wsum = wsum + wt[ordinate];
    
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wt[ordinate] = wt[ordinate]*(norm/wsum);

    // Sort the directions by xi and then by mu
    sortOrdinates();
    
    // Verify that the quadrature meets our integration requirements.
    Ensure( soft_equiv(iDomega(),norm) );

    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    Ensure( soft_equiv(iod[1],0.0) );
		    
    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( soft_equiv(iood[0],norm/3.0) );  // mu*mu
    Ensure( soft_equiv(iood[1],0.0) ); // mu*eta
    Ensure( soft_equiv(iood[2],0.0) ); // eta*mu
    Ensure( soft_equiv(iood[3],norm/3.0) ); // eta*eta

    // Copy quadrature data { mu, eta } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims = dimensionality();
    for ( size_t ordinate = 0; ordinate < numOrdinates; ++ordinate )
    {
	omega[ordinate].resize(ndims);
	omega[ordinate][0] = mu[ordinate];
	omega[ordinate][1] = xi[ordinate];
    }

    //display();

} // end of Q2DTriChebyshevLegendre constructor

//---------------------------------------------------------------------------//
/*!
 * \brief Resort all of the ordinates by xi and then by mu.
 *
 * The ctor for OrdinateSet sorts automatically.
 */
void Q2DTriChebyshevLegendre::sortOrdinates(void)
{
    size_t len( mu.size() );

    // temporary storage
    vector<Ordinate> omega;
    for( size_t m=0; m<len; ++m )
    {
        double eta=std::sqrt(1.0-mu[m]*mu[m]-xi[m]*xi[m]);

        //omega.push_back( Ordinate(mu[m],eta*(xi[m]<0?-1:1),xi[m],wt[m] ) );
        omega.push_back( Ordinate(mu[m],eta,xi[m],wt[m] ) );
    }
    
    std::sort(omega.begin(),omega.end(),Ordinate::SnCompare);
    
    // Save sorted data
    for( size_t m=0; m<len; ++m )
    {
        mu[m]=omega[m].mu();
        //xi[m]=omega[m].eta();
        xi[m]=omega[m].xi();
        wt[m]=omega[m].wt();        
    }
    
    return;
}

//---------------------------------------------------------------------------//

void Q2DTriChebyshevLegendre::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t    xi        \t     wt      " << endl;
    cout << "  --- \t------------- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ordinate = 0; ordinate < mu.size(); ++ordinate )
    {
	cout << "   "
	     << ordinate << "\t"
	     << setprecision(10) << mu[ordinate]  << "\t"
	     << setprecision(10) << xi[ordinate] << "\t"
	     << setprecision(10) << wt[ordinate]  << endl;
	sum_wt += wt[ordinate];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q2DTriChebyshevLegendre.cc
//---------------------------------------------------------------------------//
