//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/GeneralQuadrature.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  
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
#include "GeneralQuadrature.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 3D Level Symmetric quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder+2)*snOrder.
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 4*PI).
 */
GeneralQuadrature::GeneralQuadrature(size_t sn_order_,
                                     double norm_,
                                     Quadrature::QIM qm_,
                                     std::vector<double> const & mu_,
                                     std::vector<double> const & eta_,
                                     std::vector<double> const & xi_,
                                     std::vector<double> const & wt_,
                                     size_t levels_,
                                     size_t dim_,
                                     std::string const & quadratureName_,
                                     Quadrature_Class quadratureClass_) 
    : Quadrature( sn_order_, norm_, qm_ ), 
      numOrdinates ( mu_.size() ),
      numLevels ( levels_ ),
      numDims   ( dim_ ),
      quadratureName ( quadratureName_ ),
      quadratureClass ( quadratureClass_ )
{ 
    using rtt_dsxx::soft_equiv;
    using rtt_units::PI;

    Require( mu.size() == eta.size() );
    Require( mu.size() == xi.size() );
    Require( mu.size() == wt.size() );
    
    // Force the direction vectors to be the correct length.
    mu = mu_;
    eta = eta_;
    xi = xi_;
    wt = wt_;

    // Verify that the quadrature meets our integration requirements.
    Ensure( soft_equiv(iDomega(),norm) );

    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    // In 2d, eta may be a dependent variable (all entries are positive).
    // Thus, only one of the following integrals must be zero.
    if( numDims == 2 )
    {
        Ensure(soft_equiv(iod[1],0.0) || soft_equiv(iod[2],0.0) );
    }
    else
    {
        Ensure( soft_equiv(iod[1],0.0) );
        Ensure( soft_equiv(iod[2],0.0) );
    }
    
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

    // Copy quadrature data { mu, eta, xi } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims(0);
    // Set the dimensionality based on the available data.
    if(  mu.size() > 0 ) ndims++;
    if( eta.size() > 0 ) ndims++;
    if(  xi.size() > 0 ) ndims++;
    
    if( ndims == 1 )
    {
        for ( size_t i=0; i<numOrdinates; ++i )
        {
            omega[i].resize( ndims );
            omega[i][0] = mu[i];
        }
    }
    else if( ndims == 2 )
    {
        for ( size_t i=0; i<numOrdinates; ++i )
        {
            omega[i].resize( ndims );
            omega[i][0] = mu[i];
            omega[i][1] = eta[i];
        }
    }
    else if( ndims ==3 )
    {
        for ( size_t i=0; i<numOrdinates; ++i )
        {
            omega[i].resize( ndims );
            omega[i][0] = mu[i];
            omega[i][1] = eta[i];
            omega[i][2] = xi[i];
        }
    }
    else
    {
        Insist(ndims<3 && ndims>0,"Unexpected value for ndims in GeneralQuadrature.cc");
    }
    return;
} // end of Q3LevelSym() constructor.

//---------------------------------------------------------------------------//

void GeneralQuadrature::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << name() << " quadrature." 
         << endl << " Directions and weights are:" 
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
//                 end of GeneralQuadrature.cc
//---------------------------------------------------------------------------//
