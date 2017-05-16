//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/utilities.hh
 * \author Kent G. Budge
 * \brief  Declarations of a number of useful parsing utilities.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file declares functions that parse certain common constructs in a
 * uniform way.
 */
//---------------------------------------------------------------------------//

#ifndef parser_utilities_hh
#define parser_utilities_hh

#include "Expression.hh"
#include "mesh_element/Geometry.hh"
#include <memory>

namespace rtt_parser {
//! Can the next token in the stream be interpreted as real number?
DLL_PUBLIC_parser bool at_real(Token_Stream &tokens);

//! Is the next token in the stream a unit name?
DLL_PUBLIC_parser bool at_unit_term(Token_Stream &tokens,
                                    unsigned position = 0);

DLL_PUBLIC_parser unsigned parse_positive_integer(Token_Stream &);

DLL_PUBLIC_parser unsigned parse_unsigned_integer(Token_Stream &);

DLL_PUBLIC_parser int parse_integer(Token_Stream &);

DLL_PUBLIC_parser double parse_real(Token_Stream &);

DLL_PUBLIC_parser double parse_positive_real(Token_Stream &);

DLL_PUBLIC_parser double parse_nonnegative_real(Token_Stream &);

DLL_PUBLIC_parser bool parse_bool(Token_Stream &);

DLL_PUBLIC_parser Unit parse_unit(Token_Stream &);

DLL_PUBLIC_parser void parse_vector(Token_Stream &, double[]);

//! parser a quote-delimited string, stripping the quotes.
DLL_PUBLIC_parser std::string parse_manifest_string(Token_Stream &tokens);

DLL_PUBLIC_parser void
parse_geometry(Token_Stream &tokens,
               rtt_mesh_element::Geometry &parsed_geometry);

DLL_PUBLIC_parser void parse_unsigned_vector(Token_Stream &, unsigned[],
                                             unsigned);
DLL_PUBLIC_parser void
set_internal_unit_system(rtt_units::UnitSystem const &units);
DLL_PUBLIC_parser void set_unit_expressions_are_required(bool);
DLL_PUBLIC_parser rtt_units::UnitSystem const &get_internal_unit_system();
DLL_PUBLIC_parser bool unit_expressions_are_required();
DLL_PUBLIC_parser void free_internal_unit_system();

//! parser a real number followed by a unit expression.
DLL_PUBLIC_parser double parse_quantity(Token_Stream &tokens, Unit const &unit,
                                        char const *name);

//! parse an expression followed by a unit expression.
DLL_PUBLIC_parser std::shared_ptr<Expression>
parse_quantity(Token_Stream &tokens, Unit const &unit, char const *name,
               unsigned number_of_variables,
               std::map<string, pair<unsigned, Unit>> const &);

DLL_PUBLIC_parser double parse_temperature(Token_Stream &);

DLL_PUBLIC_parser std::shared_ptr<Expression>
parse_temperature(Token_Stream &, unsigned number_of_variables,
                  std::map<string, pair<unsigned, Unit>> const &);

//----------------------------------------------------------------------------//
/*! Template for parse function that produces a class object.
 *
 * The parse_class template function is intended for general use as a common
 * signature for all functions that parse an object of a specified class from a
 * Token_Stream, so that the common code idiom for parsing an object of MyClass
 * is:
 * \code
 *   std::shared_ptr<MyClass> spMyClass = parse_class<MyClass>(tokens);
 * \endcode
 * Developers may specialize this function as needed. A particular
 * implementation is suggested in Class_Parse_Table.hh.
 *
 * The template parameter names the type of the object for which a smart pointer
 * is returned.
 *
 * \param tokens Token stream from which to parse the user input.
 *
 * \return A pointer to an object matching the user specification, or NULL if
 * the specification is not valid.
 */
template <typename Class>
std::shared_ptr<Class> parse_class(Token_Stream &tokens);

//----------------------------------------------------------------------------//
/*! Template for parse function that produces a class object.
 *
 * This function resembles the previous one, but takes a second argument supply
 * a context that may alter the parser behavior..
 *
 * \param tokens Token stream from which to parse the user input.
 * \param context Context for the parse.
 * \return A pointer to an object matching the user specification, or NULL if
 * the specification is not valid.
 */
template <typename Class, typename Context>
std::shared_ptr<Class> parse_class(Token_Stream &tokens,
                                   Context const &context);

} // rtt_parser

#endif
// parser_utilities_hh

//---------------------------------------------------------------------------//
// end of utilities.hh
//---------------------------------------------------------------------------//
