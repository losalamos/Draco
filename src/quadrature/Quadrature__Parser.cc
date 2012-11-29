//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Quadrature__Parser.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Parser for various quadrature classes.
 * \note   Copyright © 2000-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "parser/Class_Parser.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "Gauss_Legendre.hh"
#include "Level_Symmetric.hh"
#include "Tri_Chebyshev_Legendre.hh"
#include "Square_Chebyshev_Legendre.hh"
#include "Double_Gauss.hh"
#include "General_Octant_Quadrature.hh"
#include "Lobatto.hh"

using namespace rtt_parser;
using namespace rtt_quadrature;

namespace // anonymous
{
//---------------------------------------------------------------------------//
SP<Quadrature> parse_gauss_legendre(Token_Stream &tokens)
{
    return Gauss_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_level_symmetric(Token_Stream &tokens)
{
    return Level_Symmetric::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_tri_cl(Token_Stream &tokens)
{
    return Tri_Chebyshev_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_square_cl(Token_Stream &tokens)
{
    return Square_Chebyshev_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_double_gauss(Token_Stream &tokens)
{
    return Double_Gauss::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_general_octant_quadrature(Token_Stream &tokens)
{
    return General_Octant_Quadrature::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parse_lobatto(Token_Stream &tokens)
{
    return Lobatto::parse(tokens);
}

//---------------------------------------------------------------------------//
SP<Quadrature> parsed_quadrature;

//---------------------------------------------------------------------------//
SP<Quadrature> &get_parsed_object()
{
    return parsed_quadrature;
}


}  // end anonymous

namespace rtt_parser
{
//---------------------------------------------------------------------------//
template<>
rtt_parser::Parse_Table
Class_Parser<Quadrature>::parse_table_(NULL, 0, Parse_Table::ONCE);

} // end namespace anonymous


namespace // anonymous
{
//---------------------------------------------------------------------------//
Parse_Table &get_parse_table()
{
    return Class_Parser<Quadrature>::parse_table_;
}

} // end anonymous

namespace rtt_parser
{
//---------------------------------------------------------------------------//
template<>
void Class_Parser<Quadrature>::post_sentinels_()
{
    // Be sure parse table has standard models
    static bool first_time = true;
    if (first_time)
    {
        Quadrature::register_quadrature("gauss legendre",
                                        parse_gauss_legendre);

        Quadrature::register_quadrature("level symmetric",
                                        parse_level_symmetric);

        Quadrature::register_quadrature("tri cl",
                                        parse_tri_cl);

        Quadrature::register_quadrature("square cl",
                                        parse_square_cl);

        Quadrature::register_quadrature("double gauss",
                                        parse_double_gauss);

        Quadrature::register_quadrature("general octant quadrature",
                                        parse_general_octant_quadrature);

        Quadrature::register_quadrature("lobatto",
                                        parse_lobatto);

        first_time = false;
    }
    
    parsed_quadrature = NULL;
}
//---------------------------------------------------------------------------//
template<>
void Class_Parser<Quadrature>::
check_completeness_(Token_Stream &tokens)
{
    tokens.check_semantics(parsed_quadrature != SP<Quadrature>(),
                           "no quadrature specified");
}

//---------------------------------------------------------------------------//
template<>
SP<Quadrature>
Class_Parser<Quadrature>::create_object_()
{
    return parsed_quadrature;
}

} //end namespace rtt_parser

namespace rtt_quadrature
{
/* static */
SP<Quadrature> Quadrature::parse(Token_Stream &tokens)
{
    Token token = tokens.shift();
    tokens.check_syntax(token.text()=="type", "expected type keyword");
    
    return Class_Parser<Quadrature>::parse(tokens);
}

//---------------------------------------------------------------------------//
/*!
 *
 * This function allows local developers to add new transport models to
 * Capsaicin, by adding their names to the set of transport models recognized
 * by the \c Transport_Model parser.
 *
 * \param keyword Name of the new transport model. Must be unique.
 *
 * \param parse_function Function to call to parse a specification for the
 * new transport model. Must return an instance of the transport model class.
 */
/* static */
void
Quadrature::register_quadrature(string const &keyword,
                                SP<Quadrature> parse_function(Token_Stream&) )
{
    Abstract_Class_Parser<Quadrature,
                          get_parse_table,
                          get_parsed_object>::register_child(keyword,
                                                             parse_function);
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------------------//
