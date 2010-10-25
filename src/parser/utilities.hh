//----------------------------------*-C++-*----------------------------------//
/*!
 * \file utilities.hh
 * \author Kent G. Budge
 * \brief Declarations of a number of useful parsing utilities.
 * \note   Copyright © 2006-2007 Los Alamos National Security, LLC
 *
 * This file declares functions that parse certain common constructs in a
 * uniform way.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Token_Stream.hh"
#include "Unit.hh"
#include "mesh_element/Geometry.hh"

namespace rtt_parser 
{
//! Can the next token in the stream be interpreted as real number?
bool at_real(Token_Stream &tokens);

unsigned parse_positive_integer(Token_Stream &);

unsigned parse_unsigned_integer(Token_Stream &);

int parse_integer(Token_Stream &);

double parse_real(Token_Stream &);

double parse_positive_real(Token_Stream &);

double parse_nonnegative_real(Token_Stream &);

bool parse_bool(Token_Stream &);

Unit parse_unit(Token_Stream &);

void parse_vector(Token_Stream &, double[]);

void parse_unsigned_vector(Token_Stream &, unsigned[], unsigned);

//! parser a real number followed by a unit expression.
double parse_quantity(Token_Stream &tokens,
		      Unit const &unit,
		      char const *name);

double parse_temperature(Token_Stream &);

//! parser a quote-delimited string, stripping the quotes.
std::string parse_manifest_string(Token_Stream &tokens);

void parse_geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry);

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of utilities.hh
//---------------------------------------------------------------------------//
