//----------------------------------*-C++-*----------------------------------//
/*!
 * \file utilities.hh
 * \author Kent G. Budge
 * \date Wed May 21 08:01:27 MDT 2003
 * \brief Declarations of a number of useful parsing utilities.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * This file declares functions that parse certain common constructs in a
 * uniform way.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Token_Stream.hh"
#include "Unit.hh"

namespace rtt_parser 
{
bool At_Real(Token_Stream &tokens);

//! Parse a positive integer.
unsigned Parse_Positive_Integer(Token_Stream &);

//! Parse an unsigned integer.
unsigned Parse_Unsigned_Integer(Token_Stream &);

//! Parse an integer.
int Parse_Integer(Token_Stream &);

double Parse_Real(Token_Stream &);

//! Parse a unit expression.
Unit Parse_Unit(Token_Stream &);

void Parse_Vector(Token_Stream &, double[]);

void Parse_Unsigned_Vector(Token_Stream &, unsigned[]);

double Parse_Quantity(Token_Stream &tokens,
		      Unit const &unit,
		      char const *name);

double Parse_Temperature(Token_Stream &);

//! Parse a quote-delimited string, stripping the quotes.
std::string Parse_Manifest_String(Token_Stream &tokens);

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of utilities.hh
//---------------------------------------------------------------------------//
