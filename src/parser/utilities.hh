//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/utilities.hh
 * \author Kent G. Budge
 * \brief  Declarations of a number of useful parsing utilities.
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file declares functions that parse certain common constructs in a
 * uniform way.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Expression.hh"
#include "Token_Stream.hh"
#include "Unit.hh"
#include "mesh_element/Geometry.hh"

namespace rtt_parser 
{
DLL_PUBLIC void disable_unit_expressions();
DLL_PUBLIC bool are_unit_expressions_disabled();

//! Can the next token in the stream be interpreted as real number?
DLL_PUBLIC bool at_real(Token_Stream &tokens);

DLL_PUBLIC unsigned parse_positive_integer(Token_Stream &);

DLL_PUBLIC unsigned parse_unsigned_integer(Token_Stream &);

DLL_PUBLIC int parse_integer(Token_Stream &);

DLL_PUBLIC double parse_real(Token_Stream &);

DLL_PUBLIC double parse_positive_real(Token_Stream &);

DLL_PUBLIC double parse_nonnegative_real(Token_Stream &);

DLL_PUBLIC bool parse_bool(Token_Stream &);

DLL_PUBLIC Unit parse_unit(Token_Stream &);

DLL_PUBLIC void parse_vector(Token_Stream &, double[]);

DLL_PUBLIC void parse_unsigned_vector(Token_Stream &, unsigned[], unsigned);

//! parser a real number followed by a unit expression.
DLL_PUBLIC double parse_quantity(
    Token_Stream &tokens,
	Unit const &unit,
	char const *name);

//! parse an expression followed by a unit expression.
DLL_PUBLIC 
SP<Expression> parse_quantity(Token_Stream &tokens,
                              Unit const &unit,
                              char const *name,
                              unsigned number_of_variables,
                              std::map<string, pair<unsigned, Unit> > const &);

DLL_PUBLIC double parse_temperature(Token_Stream &);

DLL_PUBLIC SP<Expression>
parse_temperature(Token_Stream &,
                  unsigned number_of_variables,
                  std::map<string, pair<unsigned, Unit> > const &);

//! parser a quote-delimited string, stripping the quotes.
DLL_PUBLIC std::string parse_manifest_string(Token_Stream &tokens);

DLL_PUBLIC 
void parse_geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry);

} // rtt_parser
//---------------------------------------------------------------------------//
// end of utilities.hh
//---------------------------------------------------------------------------//
