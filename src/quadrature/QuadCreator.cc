//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 15:38:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
using std::cerr;
using std::endl;

#include "Quadrature.hh"
#include "QuadCreator.hh"

#include "ds++/SP.hh"

const double PI = 3.1415926535899323846;
const double TOL = 1.0e-14;

namespace rtt_quadrature
{

// QuadCreate returns a smart pointer to a Quadrature object.  The Quadrature 
// object is created according to the clients specifications.
//
// sn_order defaults to 4 (see QuadCreator.hh)
// norm     defaults to 4*PI, 2*PI or 2 for 3D, 2D or 1D problems,
//          respectively -- See the constructors in Quadrature.cc

// if norm was not specifed by the client it will be equal to zero here and
// will be set to an appropriate default value here.

rtt_dsxx::SP<Quadrature> 
QuadCreator::QuadCreate( QuadCreator::Qid quad_type, 
			 int sn_order, double norm ) 
{
    switch( quad_type ) {
    case( GaussLeg ):
	// if the client did not specify a value for norm then it will be
	// zero here.  We must set it to a default value of 2.0.
	if ( fabs(norm) <= TOL ) norm = 2.0;
	return rtt_dsxx::SP<Quadrature> 
	    ( new Q1DGaussLeg( sn_order, norm ));
    case( LevelSym ):
	if ( fabs(norm) <= TOL ) norm = 4.0*PI;
	return rtt_dsxx::SP<Quadrature>
	    ( new Q3DLevelSym( sn_order, norm ));
    default:
	cerr << endl << "QuadCreator::QuadCreate"
	     << endl << "   --> No value for quad_type specfied.  Aborting."
	     << endl;
	return rtt_dsxx::SP<Quadrature> ();
    }
}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of Quadrature.cc
//---------------------------------------------------------------------------//
