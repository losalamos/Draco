//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadCreator.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 15:38:56 2000
 * \brief  \link rtt_quadrature::QuadCreator QuadCreator \endlink
 *         class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// Revision History:
// ------------------
// 1.4) The member function "QuadCreate" was renamed "quadCreate".
//      Moved "using" statements inside of namespace block.
//      Added "using std::fabs".
//      Explicitly create SP<Quadrature> spQuad before switch
//         statement.
//      Removed redundant declarations of PI and TOL.
//      Added Insist() call when an unknown quad_type is specified.

#include <iostream>
#include <cmath>

#include "units/PhysicalConstants.hh"

#include "Q1DGaussLeg.hh"
#include "Q2DLevelSym.hh"
#include "Q3DLevelSym.hh"
#include "QuadCreator.hh"

namespace rtt_quadrature
{

using std::cerr;
using std::endl;
using std::fabs;

// quadCreate returns a smart pointer to a Quadrature object.  The Quadrature 
// object is created according to the client's specifications.
//
// sn_order defaults to 4 (see QuadCreator.hh)
// norm     defaults ot 0.0 if the user doesn't specify a value.  If
//          norm is zero it will be set to an appropriate value before 
//          the quadrature object is created (4*pi for a 3D set, 2*pi
//          for a 2D set, or 2 for a 1D set).

rtt_dsxx::SP<Quadrature> 
QuadCreator::quadCreate( QuadCreator::Qid quad_type, 
			 size_t sn_order, double norm ) 
{
    rtt_dsxx::SP<Quadrature> spQuad;

    switch( quad_type ) 
	{
	case GaussLeg:
	    // if the client did not specify a value for norm then it will be
	    // zero here.  We must set it to a default value of 2.0.
	    if ( fabs(norm) <= TOL ) norm = 2.0;
	    spQuad = new Q1DGaussLeg( sn_order, norm );
	    break;
	    
	case LevelSym2D:
	    if ( fabs(norm) <= TOL ) norm = 2.0*rtt_units::PI;
	    spQuad = new Q2DLevelSym( sn_order, norm );
	    break;
	    
	case LevelSym:
	    if ( fabs(norm) <= TOL ) norm = 4.0*rtt_units::PI;
	    spQuad = new Q3DLevelSym( sn_order, norm );
	    break;
	    
	default:
	    Insist ( false, "Unknown value for quad_type." );
	    break;
	}
    
    return spQuad;

}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of QuadCreator.cc
//---------------------------------------------------------------------------//
