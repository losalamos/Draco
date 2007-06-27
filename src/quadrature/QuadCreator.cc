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
#include <cctype>
#include <algorithm>

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
#include "QuadServices.hh"

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
 *    type  = square CL
 *    order = 40
 *    interpolation algorithm = SVD
 *    end
 * end
 * \endcode
 * It is the responsibility of the digraph parser/creator function to identify
 * the keyword \c angle \c quadrature and make a call to this function.
 *
 * \par tokens A Token_Stream that provide text information about the quadrature set to be created.
 * \return Smart pointer to a quadrature object.
 */

struct func {
    int operator()(int x) { return std::tolower(x); }
};

rtt_dsxx::SP<Quadrature> 
QuadCreator::quadCreate( rtt_parser::Token_Stream &tokens )
{
    using namespace rtt_parser;
    using std::string;


    struct func tl;

    // Items that refine the quadrature set definition.

    QuadCreator::Qid quad_type( QuadCreator::LevelSym2D );
    double quad_norm(     4.0*rtt_units::PI );
    unsigned sn_order(    2 );     // default
    QIM      interpModel( SN );    // default

    while( tokens.Lookahead().Type() != END )
    {
        // Get next token
        Token const token = tokens.Shift();

        if (token.Type() != KEYWORD)
        {
            tokens.Report_Syntax_Error("expected a keyword");
        }

        std::string tokenText = token.Text();
        std::transform(tokenText.begin(),tokenText.end(),
                       tokenText.begin(),tl);

        
        if( tokenText == "type" )
        {
            // Get the type.
            string qtype = tokens.Shift().Text();
            // convert to use all lower case
            std::transform(qtype.begin(),qtype.end(),
                           qtype.begin(),tl);

            qidm::const_iterator pos = Qid_map.find( qtype );
            
            if( pos == Qid_map.end() )
                tokens.Report_Semantic_Error(
                    "I don't know anything about the quadrature type = "
                    +qtype);
            else
                quad_type = pos->second;
            
            // Set the default norm value
            quad_norm = norm_map.find( quad_type )->second;
        }

        // This block parses the quad_type when "type =" was not provided.
        else if( Qid_map.find( tokenText ) != Qid_map.end() )
        {
            qidm::const_iterator pos = Qid_map.find( tokenText );
            
            if( pos == Qid_map.end() )
                tokens.Report_Semantic_Error(
                    "I don't know anything about the quadrature type = "
                    +tokenText);
            else
                quad_type = pos->second;
            
            // Set the default norm value
            quad_norm = norm_map.find( quad_type )->second;
        }

        
        else if( token.Text() == "order")
        {
            sn_order = Parse_Positive_Integer(tokens);
            if (sn_order%2 != 0)
            {
                tokens.Report_Semantic_Error("quadrature order must be even");
                sn_order = 2;
            }
        }   
        else if( token.Text() == "interpolation algorithm")
        {
            string s = tokens.Shift().Text();
            std::cout << s << std::endl;
            // force lower case
            std::transform(s.begin(),s.end(),s.begin(),tl);
            if( s == "sn" )
                interpModel = SN;
            else if( s == "galerkin" )
                interpModel = GALERKIN;
            else if( s == "svd" )
                interpModel = SVD;
            else
                tokens.Report_Semantic_Error(
                    string("I don't know anything about \"angle quadrature: ")
                    +string("interpolation algorithm = ")+s
                    +string("\". Expecting one of: SN, Galerkin or SVD"));
        }
        else
        {
            tokens.Report_Syntax_Error("unrecognized keyword.  Expected \"end,\" \"order,\" or \"interpolation algorithm.\"");
        }

    } // end while

    // Read the "end" that signifies the end of the quadrature block
    Token const token = tokens.Shift();
    Ensure( token.Type() == END );
    
    rtt_dsxx::SP<Quadrature> parsed_quadrature =
        quadCreate(quad_type,sn_order, quad_norm);

    if (parsed_quadrature == rtt_dsxx::SP<rtt_quadrature::Quadrature>())
        tokens.Report_Semantic_Error("Could not construct quadrature");

    return parsed_quadrature;
}

/*!
 * \brief Generate a map that returns the appropriate enum given a string."
 */
QuadCreator::qidm QuadCreator::createQidMap(void) const
{
    qidm qid_map;

    qid_map["square cl"] = SquareCL;
    qid_map["gauss legendre"] = GaussLeg;
    qid_map["lobatto"] = Lobatto;
    qid_map["double gauss"] = DoubleGauss;
    qid_map["level symmetric"] = LevelSym2D;
    qid_map["level symmetric 3D"] = LevelSym;
    qid_map["axial"] = Axial1D;
    
    return qid_map;
}

/*!
 * \brief Generate a map that returns the appropriate norm given a quadrature type."
 */
QuadCreator::normmap QuadCreator::createNormMap(void) const
{
    normmap nmap;
    double const fourpi = 4.0*rtt_units::PI;

    nmap[SquareCL]=fourpi;
    nmap[GaussLeg]=2.0;
    nmap[Lobatto]=2.0;
    nmap[DoubleGauss]=2.0;
    nmap[LevelSym2D]=fourpi;
    nmap[LevelSym]=fourpi;
    nmap[Axial1D] = 2.0;
    
    return nmap;
}


} // end namespace rtt_quadrature


//---------------------------------------------------------------------------//
//                              end of QuadCreator.cc
//---------------------------------------------------------------------------//
