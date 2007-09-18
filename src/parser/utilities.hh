//----------------------------------*-C++-*----------------------------------//
/*!
 * \file utilities.hh
 * \author Kent G. Budge
 * \date Wed May 21 08:01:27 MDT 2003
 * \brief Declarations of a number of useful parsing utilities.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
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
bool At_Real(Token_Stream &tokens);

unsigned Parse_Positive_Integer(Token_Stream &);

unsigned Parse_Unsigned_Integer(Token_Stream &);

int Parse_Integer(Token_Stream &);

double Parse_Real(Token_Stream &);

double Parse_Positive_Real(Token_Stream &);

//! Parse a unit expression.
Unit Parse_Unit(Token_Stream &);

void Parse_Vector(Token_Stream &, double[]);

//! Parse an unsigned integer vector
void Parse_Unsigned_Vector(Token_Stream &, unsigned[], unsigned);

//! Parse a real number followed by a unit expression.
double Parse_Quantity(Token_Stream &tokens,
		      Unit const &unit,
		      char const *name);

double Parse_Temperature(Token_Stream &);

//! Parse a quote-delimited string, stripping the quotes.
std::string Parse_Manifest_String(Token_Stream &tokens);

void Parse_Geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry);

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of utilities.hh
//---------------------------------------------------------------------------//
