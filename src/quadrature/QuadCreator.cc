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

namespace rtt_quadrature
{

// QuadCreate returns a pointer to a Quadrature object.  The Quadrature
// object is created according to the clients specifications.
//
// sn_order defaults to 4 (see QuadCreator.hh)
// norm     defaults to 4*PI

Quadrature* QuadCreator::QuadCreate( QuadCreator::Qid quad_type, 
				     int sn_order, double norm ) 
{

    // verify norm > 0.0
    // verify sn_order > 1
    

    switch( quad_type ) {
    case( GaussLeg ):
	return new Q1DGaussLeg( sn_order );
    case( LevelSym ):
	return new Q3DLevelSym( sn_order, norm );
// 	cerr << endl << "QuadCreator::QuadCreate"
// 	     << endl << "   --> quad_type = LevelSym not currently available."
// 	     << endl << "       aborting" << endl;
// 	return 0;
    default:
	cerr << endl << "QuadCreator::QuadCreate"
	     << endl << "   --> No value for quad_type specfied.  Aborting."
	     << endl;
	return 0;
    }
}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of Quadrature.cc
//---------------------------------------------------------------------------//
