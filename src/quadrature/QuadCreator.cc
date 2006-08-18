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
#include "ds++/SP.hh"
#include "units/PhysicalConstants.hh"
#include "parser/utilities.hh"

#include "Q1DGaussLeg.hh"
#include "Q1DLobatto.hh"
#include "Q1DDoubleGauss.hh"
#include "Q1Axial.hh"
#include "Q2DLevelSym.hh"
#include "Q3DLevelSym.hh"
#include "Q2DSquareChebyshevLegendre.hh"
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
			 size_t sn_order,
                         double norm ) 
{
    using rtt_dsxx::soft_equiv;

    rtt_dsxx::SP<Quadrature> spQuad;

    try
    {
        
    switch( quad_type ) 
	{
	case GaussLeg:
	    // if the client did not specify a value for norm then it will be
	    // zero here.  We must set it to a default value of 2.0.
	    if ( soft_equiv(norm,0.0) ) norm = 2.0;
	    spQuad = new Q1DGaussLeg( sn_order, norm );
	    break;

	case Lobatto:
	    // if the client did not specify a value for norm then it will be
	    // zero here.  We must set it to a default value of 2.0.
	    if ( soft_equiv(norm,0.0) ) norm = 2.0;
	    spQuad = new Q1DLobatto( sn_order, norm );
	    break;

	case DoubleGauss:
	    // if the client did not specify a value for norm then it will be
	    // zero here.  We must set it to a default value of 2.0.
	    if ( soft_equiv(norm,0.0) ) norm = 2.0;
	    spQuad = new Q1DDoubleGauss( sn_order, norm );
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

	case SquareCL:
	    if ( soft_equiv(norm,0.0) ) norm = 4.0*rtt_units::PI;
	    spQuad = new Q2DSquareChebyshevLegendre( sn_order, norm );
	    break;	    

	default:
	    Insist ( false, "Unknown value for quad_type." );
	    break;
	}

    }
    catch( rtt_dsxx::assertion &error)
    {
        std::cout << "ERROR: While constructing " << quad_type << ", "
                  << error.what() << std::endl;
        throw;
    }
    
    return spQuad;

}

//---------------------------------------------------------------------------//
/*!
 * \brief quadCreate constructs a Quadrature object from a Token_Stream.
 *
 * This quadrature creator only requires 1 parameter -- a
 * rtt_parser::Token_Stream.  This function generates the appropriate
 * information for a call to the default constructor from data streams in the
 * Token_Stream.
 *
 * The Token_Stream is expected to have the following data:
 *
 * 1. A text string that describes the type of quadratre.  Valid strings are: 
 *
 *   \arg \c gauss \c legendre
 *   \arg \c level \c symmetric
 *   \arg \c square \c CL
 * 
 * 2. The quadrature order, specified by the keyword \c order followed by an
 * integer value:
 *
 *    \arg \c order \c 8
 *
 * 3. The keyword \c end to specify that we are at the end of the quadrature
 * specification block.
 * 
 * Example 1:
 * \code
 * level symmetric
 *    order 4
 * end
 * \endcode
 * When called from a solver, the quadrature block is typically initialized by
 * the keyword "angle quadrature."  In a SERRANO input deck, one would expect
 * to find the following token stream within the digraph block:
 * \code
 * angle quadrature
 *    square CL
 *       order 40
 *    end
 * end
 * \endcode
 * It is the responsibility of the digraph parser/creator function to identify
 * the keyword \c angle \c quadrature and make a call to this function.
 *
 * \par tokens A Token_Stream that provide text information about the quadrature set to be created.
 * \return Smart pointer to a quadrature object.
 */
rtt_dsxx::SP<Quadrature> 
QuadCreator::quadCreate( rtt_parser::Token_Stream &tokens )
{
    using namespace rtt_parser;
    
    QuadCreator::Qid quad_type;
    double quad_norm;

    Token const token = tokens.Shift();
    if (token.Text()=="gauss legendre")
    {
        quad_type = QuadCreator::GaussLeg;
        quad_norm = 2.0;
    }
    else if (token.Text()=="level symmetric")
    {
        quad_type = QuadCreator::LevelSym2D;
        quad_norm = 4.0*rtt_units::PI;
    }
    else if (token.Text()=="square CL")
    {
        quad_type = QuadCreator::SquareCL;
        quad_norm = 4.0*rtt_units::PI;
    }
    else
    {
        tokens.Report_Syntax_Error("expected a quadrature type specification");
    }

    if (tokens.Shift().Text() != "order")
        tokens.Report_Syntax_Error("unrecognized keyword");

    unsigned sn_order = Parse_Positive_Integer(tokens);
    if (sn_order%2 != 0)
    {
        tokens.Report_Semantic_Error("quadrature order must be even");
        sn_order = 2;
    }

    // end of quadrature type block.
    if (tokens.Shift().Type() != END)
        tokens.Report_Syntax_Error("missing 'end'");

    rtt_dsxx::SP<Quadrature> parsed_quadrature =
        quadCreate(quad_type,sn_order, quad_norm);

    if (parsed_quadrature == rtt_dsxx::SP<rtt_quadrature::Quadrature>())
        tokens.Report_Semantic_Error("Could not construct quadrature");

    return parsed_quadrature;
}

} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of QuadCreator.cc
//---------------------------------------------------------------------------//
