//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Quadrature__parser.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Parsers for various quadrature classes.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC. All rights
 *         reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "Quadrature__parser.hh"

#include "Double_Gauss.hh"
#include "Gauss_Legendre.hh"
#include "General_Octant_Quadrature.hh"
#include "Level_Symmetric.hh"
#include "Lobatto.hh"
#include "Product_Chebyshev_Legendre.hh"
#include "Square_Chebyshev_Legendre.hh"
#include "Tri_Chebyshev_Legendre.hh"
#include "parser/Abstract_Class_Parser.hh"
#include "parser/utilities.hh"

using namespace rtt_parser;
using namespace rtt_quadrature;

namespace // anonymous
{
//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_gauss_legendre(Token_Stream &tokens) {
  return Gauss_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_level_symmetric(Token_Stream &tokens) {
  return Level_Symmetric::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_tri_cl(Token_Stream &tokens) {
  return Tri_Chebyshev_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_square_cl(Token_Stream &tokens) {
  return Square_Chebyshev_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_product_cl(Token_Stream &tokens) {
  return Product_Chebyshev_Legendre::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_double_gauss(Token_Stream &tokens) {
  return Double_Gauss::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_general_octant_quadrature(Token_Stream &tokens) {
  return General_Octant_Quadrature::parse(tokens);
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> parse_lobatto(Token_Stream &tokens) {
  return Lobatto::parse(tokens);
}

} // end anonymous

namespace rtt_parser {
//---------------------------------------------------------------------------------------//
Class_Parse_Table<Quadrature> *Class_Parse_Table<Quadrature>::current_;
Parse_Table Class_Parse_Table<Quadrature>::parse_table_(NULL, 0,
                                                        Parse_Table::ONCE);
SP<Quadrature> Class_Parse_Table<Quadrature>::parsed_quadrature_;

//---------------------------------------------------------------------------------------//
Class_Parse_Table<Quadrature>::Class_Parse_Table() {
  // Be sure parse table has standard models
  static bool first_time = true;
  if (first_time) {
    register_quadrature("gauss legendre", parse_gauss_legendre);

    register_quadrature("level symmetric", parse_level_symmetric);

    register_quadrature("tri cl", parse_tri_cl);

    register_quadrature("square cl", parse_square_cl);

    register_quadrature("product cl", parse_product_cl);

    register_quadrature("double gauss", parse_double_gauss);

    register_quadrature("general octant quadrature",
                        parse_general_octant_quadrature);

    register_quadrature("lobatto", parse_lobatto);

    first_time = false;
  }

  parsed_quadrature_.reset();

  current_ = this;
}

//---------------------------------------------------------------------------------------//
void Class_Parse_Table<Quadrature>::check_completeness(Token_Stream &tokens) {
  tokens.check_semantics(parsed_quadrature_ != nullptr,
                         "no quadrature specified");
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> Class_Parse_Table<Quadrature>::create_object() {
  return parsed_quadrature_;
}

//---------------------------------------------------------------------------------------//
SP<Quadrature> &Class_Parse_Table<Quadrature>::get_parsed_object() {
  return parsed_quadrature_;
}

//---------------------------------------------------------------------------------------//
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
void Class_Parse_Table<Quadrature>::register_quadrature(
    string const &keyword, SP<Quadrature> parse_function(Token_Stream &)) {
  Abstract_Class_Parser<Quadrature, get_parse_table,
                        get_parsed_object>::register_child(keyword,
                                                           parse_function);
}

//---------------------------------------------------------------------------------------//
template <>
DLL_PUBLIC_quadrature SP<Quadrature>
parse_class<Quadrature>(Token_Stream &tokens) {
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "type", "expected type keyword");

  return parse_class_from_table<Class_Parse_Table<Quadrature>>(tokens);
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------------------//
