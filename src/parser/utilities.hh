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

#ifndef parser_utilities_hh
#define parser_utilities_hh

#include "mesh_element/Geometry.hh"
#include "units/UnitSystem.hh"
#include "Expression.hh"
#include "Token_Stream.hh"
#include "Unit.hh"

namespace rtt_parser 
{
//! Can the next token in the stream be interpreted as real number?
DLL_PUBLIC bool at_real(Token_Stream &tokens);

//! Is the next token in the stream a unit name?
DLL_PUBLIC bool at_unit_term(Token_Stream &tokens, unsigned position=0);

DLL_PUBLIC unsigned parse_positive_integer(Token_Stream &);

DLL_PUBLIC unsigned parse_unsigned_integer(Token_Stream &);

DLL_PUBLIC int parse_integer(Token_Stream &);

DLL_PUBLIC double parse_real(Token_Stream &);

DLL_PUBLIC double parse_positive_real(Token_Stream &);

DLL_PUBLIC double parse_nonnegative_real(Token_Stream &);

DLL_PUBLIC bool parse_bool(Token_Stream &);

DLL_PUBLIC Unit parse_unit(Token_Stream &);

DLL_PUBLIC void parse_vector(Token_Stream &, double[]);

//! parser a quote-delimited string, stripping the quotes.
DLL_PUBLIC std::string parse_manifest_string(Token_Stream &tokens);

DLL_PUBLIC 
void parse_geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry);


DLL_PUBLIC void parse_unsigned_vector(Token_Stream &, unsigned[], unsigned);
DLL_PUBLIC void set_internal_unit_system(rtt_units::UnitSystem const &units);
DLL_PUBLIC void set_unit_expressions_are_required(bool);
DLL_PUBLIC rtt_units::UnitSystem const &get_internal_unit_system();
DLL_PUBLIC bool unit_expressions_are_required();

//! parser a real number followed by a unit expression.
DLL_PUBLIC double parse_quantity(Token_Stream &tokens,
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

} // rtt_parser

#endif
// parser_utilities_hh

//---------------------------------------------------------------------------//
// end of utilities.hh
//---------------------------------------------------------------------------//
