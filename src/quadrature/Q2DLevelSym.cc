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
#include <algorithm>

#include "ds++/Soft_Equivalence.hh"
#include "Q3DLevelSym.hh"
#include "Q2DLevelSym.hh"
#include "Ordinate.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructs a 2D Level Symmetric quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder+2)*snOrder/2.
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2*PI).
 */
//---------------------------------------------------------------------------//
Q2DLevelSym::Q2DLevelSym( size_t sn_order_, double norm_ ) 
    : Quadrature( sn_order_, norm_ ), numOrdinates (sn_order_ * (sn_order_+2)/2)
{
    using std::fabs;
    using rtt_dsxx::soft_equiv;

    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    Require ( snOrder%2 == 0 );

    if ( snOrder > 24 )
    {
        throw std::range_error("Sn order out of range");
    }
    
    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates);
    xi.resize(numOrdinates);
    wt.resize(numOrdinates);

    // Use the 3D level-symmetric quadrature to construct the 2D one.
    // In the future we probably should reverse the dependence - we should
    // use the 2D quadrature to construct the 3D one, rather than depending
    // on a particular data layout in the 3D quadrature.
    Q3DLevelSym quad3D(snOrder, 2.0*norm);
    for( size_t ordinate = 0; ordinate < numOrdinates; ++ordinate )
    {
	mu[ordinate] = quad3D.getMu(ordinate);
	xi[ordinate] = quad3D.getEta(ordinate);
	wt[ordinate] = quad3D.getWt(ordinate);
    }

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

#ifdef ENSURE_ON
    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    Ensure( soft_equiv(iod[1],0.0) );
		    
    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( soft_equiv(iood[0],norm/3.0) );  // mu*mu
    Ensure( soft_equiv(iood[1],0.0) ); // mu*xi
    Ensure( soft_equiv(iood[2],0.0) ); // xi*mu
    Ensure( soft_equiv(iood[3],norm/3.0) ); // xi*xi
#endif
    // ENSURE_ON

    // Copy quadrature data { mu, xi } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims = dimensionality();
    for ( size_t ordinate = 0; ordinate < numOrdinates; ++ordinate )
    {
	omega[ordinate].resize(ndims);
	omega[ordinate][0] = mu[ordinate];
	omega[ordinate][1] = xi[ordinate];
    }
} // end of Q2DLevelSym() constructor.

//---------------------------------------------------------------------------//
/*!
 * \brief Resort all of the ordinates by xi and then by mu.
 *
 * The ctor for OrdinateSet sorts automatically.
 */
void Q2DLevelSym::sortOrdinates(void)
{
    size_t len( mu.size() );

    // temporary storage
    vector<Ordinate> omega;
    for( size_t m=0; m<len; ++m )
    {
        double eta=std::sqrt(1.0-mu[m]*mu[m]-xi[m]*xi[m]);
        omega.push_back( Ordinate(mu[m],eta,xi[m],wt[m] ) );    
    }
    
    std::sort(omega.begin(),omega.end(),Ordinate::SnCompare);
    
    // Save sorted data
    for( size_t m=0; m<len; ++m )
    {
        mu[m]=omega[m].mu();
        xi[m]=omega[m].xi();
        wt[m]=omega[m].wt();        
    }
    
    return;
}

//---------------------------------------------------------------------------//

void Q2DLevelSym::display() const 
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
//                 end of Q2DLevelSym.cc
//---------------------------------------------------------------------------//
