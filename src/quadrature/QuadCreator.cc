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

#include <iostream>
#include <cmath>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "Q1DGaussLeg.hh"
#include "Q1Axial.hh"
#include "Q2DLevelSym.hh"
#include "Q3DLevelSym.hh"
#include "QuadCreator.hh"

namespace rtt_quadrature
{

/*!
 * \brief quadCreate constructs a Quadrature object.
 *
 * The Quad creator only requires 1 parameter -- the quadrature
 * identifier (see QuadCreator::Qid).  The two addtional parameters can
 * optionally be used to specify the sn_order and a normalization for the 
 * quadrature weights.  The sn_order defaults to 4 and the default value
 * for normalization constant varies with the dimensionality of the
 * quadrature set (2, 2*pi or 4*pi for 1D, 2D or 3D sets).
 *
 * Another parameter may need to be added to this constructor to specify
 * the number of dimensions requested.  Currently Qid directly specifies
 * the dimensionality of the quadrature set.
 *
 * \par quad_type An identifier that specifies the type of quadrature
 *                  to construct.
 * \par sn_order  The SN order for the constructed quadrature
 *                  set. (Default: 4)
 * \par norm      The sum of the quadrature weights are forced to sum
 *                  to this value. (Default: 2, 2*pi or 4*pi based on the 
 *                  dimensionality of the quadrature set.)
 * \return Smart pointer to a quadrature object.
 */
rtt_dsxx::SP<Quadrature> 
QuadCreator::quadCreate( QuadCreator::Qid quad_type, 
			 size_t sn_order, double norm ) 
{
    using rtt_dsxx::soft_equiv;

    rtt_dsxx::SP<Quadrature> spQuad;

    switch( quad_type ) 
	{
	case GaussLeg:
	    // if the client did not specify a value for norm then it will be
	    // zero here.  We must set it to a default value of 2.0.
	    if ( soft_equiv(norm,0.0) ) norm = 2.0;
	    spQuad = new Q1DGaussLeg( sn_order, norm );
	    break;
	    
	case Axial1D:
	    if ( soft_equiv(norm,0.0) ) norm = 2.0;
	    spQuad = new Q1Axial( sn_order, norm );
	    break;
	    
	case LevelSym2D:
	    if ( soft_equiv(norm,0.0) ) norm = 2.0*rtt_units::PI;
	    spQuad = new Q2DLevelSym( sn_order, norm );
	    break;
	    
	case LevelSym:
	    if ( soft_equiv(norm,0.0) ) norm = 4.0*rtt_units::PI;
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
