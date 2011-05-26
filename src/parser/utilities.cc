//----------------------------------*-C++-*----------------------------------//
/*!
 * \file utilities.cc
 * \author Kent G. Budge
 * \date 18 Feb 2003
 * \brief Definitions of parsing utility functions.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include "utilities.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_parser 
{
using namespace std;

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

unsigned parse_unsigned_integer(Token_Stream &tokens)
{
    Token const token = tokens.shift();
    if (token.type() == INTEGER)
    {
	errno = 0;
	char *endptr;
	unsigned long const Result = strtoul(token.text().c_str(), &endptr, 0);
	if (Result != static_cast<unsigned>(Result) || errno==ERANGE)
	{
	    tokens.report_semantic_error("integer value overflows");
	}
	Check(endptr != NULL);
	return Result;
    }
    else
    {
	tokens.report_syntax_error(token, "expected an unsigned integer");
	return 0;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

unsigned parse_positive_integer(Token_Stream &tokens)
{
    unsigned  Result = parse_unsigned_integer(tokens);
    if (Result==0)
    {
	tokens.report_semantic_error("expected a positive integer");
	Result = 1;
    }

    Ensure(Result>0);
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

int parse_integer(Token_Stream &tokens)
{
    Token token = tokens.shift();
    string text;
    if (token.text() == "+")
    {
        token = tokens.shift();
    }
    else if (token.text() == "-")
    {
        text = '-';
        token = tokens.shift();
    }
    if (token.type() == INTEGER)
    {
        text += token.text();
	errno = 0;
	char *endptr;
	const long Result = strtol(text.c_str(), &endptr, 0);
	if (Result != static_cast<int>(Result) || errno==ERANGE)
	{
	    tokens.report_semantic_error("integer value overflows");
	}
	Check(endptr != NULL && *endptr=='\0');
	return Result;
    }
    else
    {
	tokens.report_syntax_error(token, "expected an integer");
	return 0;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * This function does not move the cursor in the token stream.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \return \c true if the next token is REAL or INTEGER; \c false otherwise.
 */

bool at_real(Token_Stream &tokens)
{
    Token token = tokens.lookahead();
    if (token.text() == "-" || token.text() == "+")
    {
        token = tokens.lookahead(1);
    }
    return (token.type()==REAL || token.type()==INTEGER);
}

//---------------------------------------------------------------------------//
/*! 
 * We permit an integer token to appear where a real is expected, consistent
 * with the integers being a subset of reals, and with about five decades of
 * common practice in the computer language community.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */

double parse_real(Token_Stream &tokens)
{
    Token token = tokens.shift();
    string text;
    if (token.text() == "+")
    {
        token = tokens.shift();
    }
    else if (token.text() == "-")
    {
        text = '-';
        token = tokens.shift();
    }
    if (token.type() == REAL || token.type() == INTEGER)
    {
        text += token.text();
	errno = 0;
	char *endptr;
	const double Result = strtod(text.c_str(), &endptr);
	if (errno==ERANGE)
	{
	    tokens.report_semantic_error("real value overflows");
	}
	Check(endptr != NULL && *endptr=='\0');
	return Result;
    }
    else
    {
	tokens.report_syntax_error(token, "expected a real number");
	return 0.0;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

double parse_positive_real(Token_Stream &tokens)
{
    double Result = parse_real(tokens);
    if (Result<=0.0)
    {
	tokens.report_semantic_error("expected a positive quantity");
	Result = 1;
    }

    Ensure(Result>0);
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

double parse_nonnegative_real(Token_Stream &tokens)
{
    double Result = parse_real(tokens);
    if (Result<0.0)
    {
	tokens.report_semantic_error("expected a nonnegative quantity");
	Result = 1;
    }

    Ensure(Result>=0);
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param x
 * On return, contains the parsed vector components.
 * \pre \c x!=NULL
 */

void parse_vector(Token_Stream &tokens, double x[])
{
    Require(x!=NULL);

    // At least one component must be present.
    x[0] = parse_real(tokens);

    if (at_real(tokens))
    {
	x[1] = parse_real(tokens);
	if (at_real(tokens))
	{
	    x[2] = parse_real(tokens);
	}
	else
	{
	    x[2] = 0.0;
	}
    }
    else
    {
	x[1] = 0.0;
	x[2] = 0.0;
    }
}


//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param x
 * On return, contains the parsed vector components.
 * \pre \c x!=NULL
 */

void parse_unsigned_vector(Token_Stream &tokens, unsigned x[], unsigned size)
{
    Require(x!=NULL);

    for ( unsigned i = 0; i < size; ++i )
    {
        if (at_real(tokens))
        {
            x[i] = parse_unsigned_integer(tokens);
        }
        else
        {
            ostringstream str;
            str << "unsigned integer sequence too short; expected "
                << size << " values.";

            tokens.report_semantic_error(str.str().c_str());
            return;
        }
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Are we at a unit term?
 * 
 * We are at a unit term if the next token on the Token_Stream is a valid
 * unit name. 
 * 
 * \param tokens
 * Token_Stream from which to parse.
 * \param pos
 * Position in Token_Stream at which to parse.  This lookahead capability is
 * needed by parse_unit to see if a hyphen '-' is part of a unit expression.
 *
 * \return \c true if we are at the start of a unit term; \c false otherwise
 */

bool at_unit_term(Token_Stream &tokens, unsigned position = 0)
{
    Token const token = tokens.lookahead(position);
    if (token.type()==KEYWORD)
    {
	string const u = token.text();
	switch( u[0] )
	{
	case 'A':
	    // return (u[1]=='\0');
	    return (u.size() == 1);

	case 'C':
	    return (u.size() == 1);

	case 'F':
	    return (u.size() == 1);

	case 'H':
	    return (u.size() == 1 || token.text()=="Hz");

	case 'J':
	    return (u.size() == 1);

	case 'K':
	    return (u.size() == 1);

	case 'N':
	    return (u.size() == 1);

	case 'P':
	    return (token.text()=="Pa");

	case 'S':
	    return (u.size() == 1);

	case 'T':
	    return (u.size() == 1);

	case 'V':
	    return (u.size() == 1);

	case 'W':
	    return (u.size() == 1 || token.text()=="Wb");

	case 'c':
	    return (token.text()=="cd" || token.text()=="cm");

	case 'd':
	    return (token.text()=="dyne");

	case 'e':
	    return (token.text()=="erg");

	case 'f':
	    return (token.text()=="foot");

	case 'g':
	    return (u.size() == 1);

	case 'i':
	    return (token.text()=="inch");

	case 'k':
	    return (token.text()=="kg" || token.text()=="keV");

	case 'l':
	    return (token.text()=="lm" || token.text()=="lx");

	case 'm':
	    return (u.size() == 1 || token.text() == "mol");

	case 'o':
	    return (token.text()=="ohm");

	case 'p':
	    return (token.text()=="pound");

	case 'r':
	    return (token.text()=="rad");

	case 's':
	    return (u.size() == 1 || token.text()=="sr");

	default:
	    return false;
	}
    }
    else if (token.type()==OTHER && token.text()=="(")
    {
	return true;
    }
    else
    {
	return false;
    }
}

Unit parse_unit(Token_Stream &tokens);

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a unit name.
 * 
 * A unit name can either be a literal unit name, like "kg", or a
 * parenthesized unit expression.
 *
 * \param tokens
 * Token_Stream from which to parse.
 * \return The unit.
 */

static Unit parse_unit_name(Token_Stream &tokens)
{
    Token token = tokens.shift();
    if (token.type()==KEYWORD)
    {
	string const u = token.text();
	switch ( u[0] )
	{
	case 'A':
	    if ( u.size() == 1 )
		return A;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'C':
	    if ( u.size() == 1 )
		return C;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'F':
	    if ( u.size() == 1 )
		return F;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'H':
	    if ( u.size() == 1 )
		return H;
	    else if ( u.size() == 2 )
		return Hz;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'J':
	    if ( u.size() == 1 )
		return J;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'K':
	    if ( u.size() == 1 )
		return K;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'N':
	    if ( u.size() == 1 )
		return N;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'P':
	    if ( u.size() == 2 )
		return Pa;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'S':
	    if ( u.size() == 1 )
		return S;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'T':
	    if ( u.size() == 1 )
		return T;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'V':
	    if ( u.size() == 1 )
		return V;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'W':
	    if ( u.size() == 1 )
		return W;
	    else if (token.text()=="Wb")
		return Wb;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'c':
	    if (token.text()=="cd")
		return cd;
	    else if (token.text()=="cm")
		return cm;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'd':
	    if (token.text()=="dyne")
		return dyne;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'e':
	    if (token.text()=="erg")
		return erg;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'f':
	    if (token.text()=="foot")
		return foot;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'g':
	    if ( u.size() == 1 )
		return g;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'i':
	    if (token.text()=="inch")
		return inch;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'k':
	    if (token.text()=="kg")
		return kg;
	    else if (token.text()=="keV")
		return keV;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'l':
	    if (token.text()=="lm")
		return lm;
	    else if (token.text()=="lx")
		return lx;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'm':
	    if ( u.size() == 1 )
		return m;
	    else if (token.text() == "mol")
		return mol;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'o':
	    if (token.text()=="ohm")
		return ohm;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'p':
	    if (token.text()=="pound")
		return pound;
	    else
		tokens.report_syntax_error("expected a unit");

	case 'r':
	    if (token.text()=="rad")
		return rad;
	    else
		tokens.report_syntax_error("expected a unit");

	case 's':
	    if ( u.size() == 1 )
		return s;
	    else if (token.text()=="sr")
		return sr;
	    else
		tokens.report_syntax_error("expected a unit");

	default:
	    tokens.report_syntax_error("expected a unit");
	}
    }
    else if (token.type()==OTHER && token.text()=="(")
    {
	Unit Result = parse_unit(tokens);
        token = tokens.shift();
	if (token.type()!=OTHER || token.text()!=")")
	    tokens.report_syntax_error("missing ')'");
	return Result;
    }
    else
    {
	tokens.report_syntax_error("expected a unit expression");
    }
    // never reached but causes warnings
    return dimensionless;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a unit term.
 * 
 * A unit term is a unit name optionally raised to some power.
 *
 * \param tokens
 * Token_Stream from which to parse.
 * \return The unit term.
 */

static Unit parse_unit_term(Token_Stream &tokens)
{
    Unit const Result = parse_unit_name(tokens);
    Token const token = tokens.lookahead();
    if (token.text()=="^")
    {
	tokens.shift();
	double const exponent = parse_real(tokens);
	return pow(Result, exponent);
    }
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * A unit expression is a sequence of tokens with a form such as "kg-m/sec" or
 * "erg/cm^2/sec/Hz" that gives the dimensions of a physical quantity.  This
 * function parses such an expression from its Token_Stream, returning the
 * result as a Unit whose conversion factor is relative to SI.  An empty unit
 * expression is permitted and returns rtt_parser::dimensionless, the
 * identity Unit representing the pure number 1.
 * 
 * \param tokens
 * Token_Stream from which to parse.
 * \return The unit expression.
 */

Unit parse_unit(Token_Stream &tokens)
{
    if (!at_unit_term(tokens)) return dimensionless;

    Unit Result = parse_unit_term(tokens);

    for (;;)
    {
	Token const token = tokens.lookahead();
	if (token.type()==OTHER)
	{
	    if (token.text()=="-" && at_unit_term(tokens, 1))
	    {
		tokens.shift();
		Result = Result * parse_unit_term(tokens);
	    }
	    else if (token.text()=="/")
	    {
		tokens.shift();
		Result = Result / parse_unit_term(tokens);
	    }
	    else
	    {
		return Result;
	    }
	}
	else
	{
	    return Result;
	}
    }
}

//---------------------------------------------------------------------------//
/*! 
 * This function parses a quantity having dimensions. It is assumed that the
 * client expects certain dimensions for the quantity, and an exception is
 * thrown if the dimensions are not what the client expected. The quantity
 * will be converted to the desired unit system, as indicated by the \c .conv
 * member of the \c target_unit argument.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param target_unit
 * Expected units for the quantity parsed, including conversion factor.
 * \param name
 * Name of the units expected for the quantity parsed, such as "length" or
 * "ergs/cm/sec/Hz". Used to generate diagnostic messages.
 *
 * \return The parsed value, converted to the desired unit system. 
 */

double parse_quantity(Token_Stream &tokens,
		      Unit const &target_unit,
		      char const *const name)
{
    double const value = parse_real(tokens);
    Unit const unit = parse_unit(tokens);
    if (!is_compatible(unit, target_unit))
    {
	ostringstream buffer;
	buffer << "expected quantity with dimensions of " << name;
	tokens.report_semantic_error(buffer.str().c_str());
    }
    return  value*unit.conv/target_unit.conv;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a temperature specification 
 *
 * It is very common for transport researchers to specify a temperature in
 * units of energy, using Boltzmann's constant as the conversion factor.
 * This function is useful for parsers that accomodate this convention.
 * 
 * \param tokens
 * Token stream from which to parse the specification.
 *
 * \return The parsed temperature.
 *
 * \post \c Result>=0.0
 */

double parse_temperature(Token_Stream &tokens)
{
    double const T = parse_real(tokens);
    Unit const u = parse_unit(tokens);
    double Result;
    if (is_compatible(u, K))
    {
	Result = T * u.conv;
    }
    else if (is_compatible(u, J))
    {
	Result = T * u.conv/rtt_units::boltzmannSI;
    }
    else
    {
	tokens.report_syntax_error("expected quantity with units of "
				   "temperature");
	return 0.0;
    }
    if (Result<0.0)
    {
        tokens.report_semantic_error("temperature must be nonnegative");
        return 0.0;
    }
    else
    {
        return Result;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a temperature specification 
 *
 * It is very common for transport researchers to specify a temperature in
 * units of energy, using Boltzmann's constant as the conversion factor.
 * This function is useful for parsers that accomodate this convention.
 * 
 * \param tokens
 * Token stream from which to parse the specification.
 *
 * \return The parsed temperature.
 *
 * \post \c Result>=0.0
 */

SP<Expression>
parse_temperature(Token_Stream &tokens,
                  unsigned const number_of_variables,
                  std::map<string, pair<unsigned, Unit> > const &variable_map)
{
    SP<Expression> T = Expression::parse(number_of_variables,
                                         variable_map,
                                         tokens);
    
    Unit const u = parse_unit(tokens)*T->units();
    if (is_compatible(u, K))
    {
        T->set_units(u);
    }
    else if (is_compatible(u, J))
    {
        T->set_units(u*K/(J*rtt_units::boltzmannSI));
    }
    else
    {
	tokens.report_syntax_error("expected quantity with units of "
				   "temperature");
    }
    return T;
}

//---------------------------------------------------------------------------//
/*! 
 * Parses a STRING token and strips the delimiting quotation marks.
 *
 * \param tokens
 * Token_Stream from which to parse.
 * \return The stripped string.
 */

std::string parse_manifest_string(Token_Stream &tokens)
{
    Token const token = tokens.shift();
    if (token.type() != STRING)
    {
	tokens.report_syntax_error("expected a string, but saw " + 
				   token.text());
    }
    string Result = token.text();
    string::size_type const length = Result.size();
    Check(length>1);
    Result = Result.substr(1, length-2);
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a geometry specification.
 * 
 * \param tokens
 * Token stream from which to parse the geometry.
 * \param parsed_geometry On entry, if the value is not \c END_GEOMETRY, a
 * diagnostic is generated to the token stream. On return, contains the
 * geometry that was parsed.
 *
 * \post <code> parsed_geometry == rtt_mesh_element::AXISYMMETRIC ||
 *         parsed_geometry == rtt_mesh_element::CARTESIAN    ||
 *         parsed_geometry == rtt_mesh_element::SPHERICAL </code>
 */

void parse_geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry)
{
    if (parsed_geometry != rtt_mesh_element::END_GEOMETRY)
    {
        tokens.report_semantic_error("geometry specified twice");
    }
    Token const token = tokens.shift();
    if (token.text() == "axisymmetric" ||
        token.text() == "cylindrical")
    {
        parsed_geometry = rtt_mesh_element::AXISYMMETRIC;
    }
    else if (token.text() == "cartesian" ||
             token.text() == "xy" ||
             token.text() == "slab")
    {
        parsed_geometry = rtt_mesh_element::CARTESIAN;
    }
    else if (token.text() == "spherical")
	{
	    parsed_geometry = rtt_mesh_element::SPHERICAL;
	}
    else
    {
        tokens.report_syntax_error(token,
                                   "expected a geometry option, but saw " + 
                                   token.text());
    }
    Ensure(parsed_geometry == rtt_mesh_element::AXISYMMETRIC ||
           parsed_geometry == rtt_mesh_element::CARTESIAN    ||
           parsed_geometry == rtt_mesh_element::SPHERICAL);
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \param tokens
 * Token stream from which to parse the quantity.
 *
 * \return The parsed quantity.
 */

bool parse_bool(Token_Stream &tokens)
{
    Token const token = tokens.shift();
    if (token.text()=="true")
    {
        return true;
    }
    else if (token.text()=="false")
    {
        return false;
    }
    else
    {
        tokens.report_syntax_error("expected 'true' or 'false'");
        return true; // to turn off warning; never reached
    }
}

//---------------------------------------------------------------------------//
/*! 
 * This function parses an expression having dimensions. It is assumed that the
 * client expects certain dimensions for the quantity, and an exception is
 * thrown if the dimensions are not what the client expected. The quantity
 * will be converted to the desired unit system, as indicated by the \c .conv
 * member of the \c target_unit argument.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param target_unit
 * Expected units for the quantity parsed, including conversion factor.
 * \param name
 * Name of the units expected for the quantity parsed, such as "length" or
 * "ergs/cm/sec/Hz". Used to generate diagnostic messages.
 *
 * \return The parsed value, converted to the desired unit system. 
 */

SP<Expression> parse_quantity(Token_Stream &tokens,
                              Unit const &target_unit,
                              char const *const name,
                              unsigned const number_of_variables,
                              std::map<string, pair<unsigned, Unit> >
                                const &variable_map)
{
    SP<Expression> value = Expression::parse(number_of_variables,
                                             variable_map,
                                             tokens);
    
    Unit unit = parse_unit(tokens)*value->units();
    if (!is_compatible(unit, target_unit))
    {
	ostringstream buffer;
	buffer << "expected quantity with dimensions of " << name;
	tokens.report_semantic_error(buffer.str().c_str());
    }
    value->set_units(unit);
    return value;
}

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of utilities.cc
//---------------------------------------------------------------------------//
