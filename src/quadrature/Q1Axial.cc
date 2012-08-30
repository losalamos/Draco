//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q1Axial.cc
 * \author Jae H Chang
 * \date   Thu Oct  7 16:05:32 2004
 * \brief  1D Axial Quadrature - used in SF TSA
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>

#include "units/PhysicalConstants.hh"
#include "Quadrature.hh"
#include "Q1Axial.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 1D Axial quadrature object.  Used for 
 *        Axial sweeps in filter sweeps.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed. For Axial it is a constant (2).
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2.0).
 */

Q1Axial::Q1Axial( size_t n, double norm_, Quadrature::QIM qm_ )
    : Quadrature( n, norm_, qm_ ), numOrdinates( n )
{
    // We require the sn_order to be 2.
    Require( n == 2.0 );
  
    // We require the normalization constant to be 2.
    Require( norm == 2.0 );

    // size the member data vectors
    mu.resize(n);
    wt.resize(n);

    double mu1 = -1.0;
    double mu2 =  1.0;
    mu[0] = mu1;
    mu[1] = mu2;
    wt[0] = 1.0;
    wt[1] = 1.0;

} // end of Q1Axial() constructor.

//---------------------------------------------------------------------------//

void Q1Axial::display() const 
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
	     << setprecision(5)  << ix     << "\t\t"
	     << setprecision(10) << mu[ix] << "\t\t"
	     << setprecision(10) << wt[ix] << endl;
	sum_wt += wt[ix];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}
 

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q1Axial.cc
//---------------------------------------------------------------------------//
