//----------------------------------*-C++-*----------------------------------//
/*!
 * \file utilities.cc
 * \author Kent G. Budge
 * \date 18 Feb 2003
 * \brief Definitions of parsing utility functions.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
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

unsigned Parse_Unsigned_Integer(Token_Stream &tokens)
{
    Token const token = tokens.Shift();
    if (token.Type() == INTEGER)
    {
	errno = 0;
	char *endptr;
	unsigned long const Result = strtoul(token.Text().c_str(), &endptr, 0);
	if (Result != static_cast<unsigned>(Result) || errno==ERANGE)
	{
	    tokens.Report_Semantic_Error("integer value overflows");
	}
	Check(endptr != NULL && *endptr=='\0');
	return Result;
    }
    else
    {
	tokens.Report_Syntax_Error(token, "expected an unsigned integer");
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

unsigned Parse_Positive_Integer(Token_Stream &tokens)
{
    unsigned  Result = Parse_Unsigned_Integer(tokens);
    if (Result==0)
    {
	tokens.Report_Semantic_Error("expected a positive integer");
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

int Parse_Integer(Token_Stream &tokens)
{
    Token token = tokens.Shift();
    string text;
    if (token.Text() == "+")
    {
        token = tokens.Shift();
    }
    else if (token.Text() == "-")
    {
        text = '-';
        token = tokens.Shift();
    }
    if (token.Type() == INTEGER)
    {
        text += token.Text();
	errno = 0;
	char *endptr;
	const long Result = strtol(text.c_str(), &endptr, 0);
	if (Result != static_cast<int>(Result) || errno==ERANGE)
	{
	    tokens.Report_Semantic_Error("integer value overflows");
	}
	Check(endptr != NULL && *endptr=='\0');
	return Result;
    }
    else
    {
	tokens.Report_Syntax_Error(token, "expected an integer");
	return 0;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Is the next token compatible with Parse_Real?
 *
 * This function does not change the token stream.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \return \c true if the next token is REAL or INTEGER; \c false otherwise.
 */

bool At_Real(Token_Stream &tokens)
{
    Token token = tokens.Lookahead();
    if (token.Text() == "-" || token.Text() == "+")
    {
        token = tokens.Lookahead(1);
    }
    return (token.Type()==REAL || token.Type()==INTEGER);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a real quantity.
 * 
 * We permit an integer token to appear where a real is expected, consistent
 * with the integers being a subset of reals, and with about five decades of
 * common practice in the computer language community.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \return The parsed quantity.
 */

double Parse_Real(Token_Stream &tokens)
{
    Token token = tokens.Shift();
    string text;
    if (token.Text() == "+")
    {
        token = tokens.Shift();
    }
    else if (token.Text() == "-")
    {
        text = '-';
        token = tokens.Shift();
    }
    if (token.Type() == REAL || token.Type() == INTEGER)
    {
        text += token.Text();
	errno = 0;
	char *endptr;
	const double Result = strtod(text.c_str(), &endptr);
	if (errno==ERANGE)
	{
	    tokens.Report_Semantic_Error("real value overflows");
	}
	Check(endptr != NULL && *endptr=='\0');
	return Result;
    }
    else
    {
	tokens.Report_Syntax_Error(token, "expected a real number");
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

double Parse_Positive_Real(Token_Stream &tokens)
{
    double Result = Parse_Real(tokens);
    if (Result<=0.0)
    {
	tokens.Report_Semantic_Error("expected a positive quantity");
	Result = 1;
    }

    Ensure(Result>0);
    return Result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Parse a vector.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param x
 * On return, contains the parsed vector components.
 * \pre \c x!=NULL
 */

void Parse_Vector(Token_Stream &tokens, double x[])
{
    Require(x!=NULL);

    // At least one component must be present.
    x[0] = Parse_Real(tokens);

    if (At_Real(tokens))
    {
	x[1] = Parse_Real(tokens);
	if (At_Real(tokens))
	{
	    x[2] = Parse_Real(tokens);
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
 * \brief Parse a vector of unsigned.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param x
 * On return, contains the parsed vector components.
 * \pre \c x!=NULL
 */

void Parse_Unsigned_Vector(Token_Stream &tokens, unsigned x[], unsigned size)
{
    Require(x!=NULL);

    for ( int i = 0; i < size; ++i )
    {
        if (At_Real(tokens))
            x[i] = Parse_Unsigned_Integer(tokens);
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
 * needed by Parse_Unit to see if a hyphen '-' is part of a unit expression.
 *
 * \return \c true if we are at the start of a unit term; \c false otherwise
 */

bool At_Unit_Term(Token_Stream &tokens, unsigned position = 0)
{
    Token const token = tokens.Lookahead(position);
    if (token.Type()==KEYWORD)
    {
	string const u = token.Text();
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
	    return (u.size() == 1 || token.Text()=="Hz");

	case 'J':
	    return (u.size() == 1);

	case 'K':
	    return (u.size() == 1);

	case 'N':
	    return (u.size() == 1);

	case 'P':
	    return (token.Text()=="Pa");

	case 'S':
	    return (u.size() == 1);

	case 'T':
	    return (u.size() == 1);

	case 'V':
	    return (u.size() == 1);

	case 'W':
	    return (u.size() == 1 || token.Text()=="Wb");

	case 'c':
	    return (token.Text()=="cd" || token.Text()=="cm");

	case 'd':
	    return (token.Text()=="dyne");

	case 'e':
	    return (token.Text()=="erg");

	case 'f':
	    return (token.Text()=="foot");

	case 'g':
	    return (u.size() == 1);

	case 'i':
	    return (token.Text()=="inch");

	case 'k':
	    return (token.Text()=="kg" || token.Text()=="keV");

	case 'l':
	    return (token.Text()=="lm" || token.Text()=="lx");

	case 'm':
	    return (u.size() == 1 || token.Text() == "mol");

	case 'o':
	    return (token.Text()=="ohm");

	case 'p':
	    return (token.Text()=="pound");

	case 'r':
	    return (token.Text()=="rad");

	case 's':
	    return (u.size() == 1 || token.Text()=="sr");

	default:
	    return false;
	}
    }
    else if (token.Type()==OTHER && token.Text()=="(")
    {
	return true;
    }
    else
    {
	return false;
    }
}

Unit Parse_Unit(Token_Stream &tokens);

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

static Unit Parse_Unit_Name(Token_Stream &tokens)
{
    Token token = tokens.Shift();
    if (token.Type()==KEYWORD)
    {
	string const u = token.Text();
	switch ( u[0] )
	{
	case 'A':
	    if ( u.size() == 1 )
		return A;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'C':
	    if ( u.size() == 1 )
		return C;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'F':
	    if ( u.size() == 1 )
		return F;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'H':
	    if ( u.size() == 1 )
		return H;
	    else if ( u.size() == 2 )
		return Hz;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'J':
	    if ( u.size() == 1 )
		return J;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'K':
	    if ( u.size() == 1 )
		return K;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'N':
	    if ( u.size() == 1 )
		return N;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'P':
	    if ( u.size() == 2 )
		return Pa;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'S':
	    if ( u.size() == 1 )
		return S;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'T':
	    if ( u.size() == 1 )
		return T;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'V':
	    if ( u.size() == 1 )
		return V;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'W':
	    if ( u.size() == 1 )
		return W;
	    else if (token.Text()=="Wb")
		return Wb;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'c':
	    if (token.Text()=="cd")
		return cd;
	    else if (token.Text()=="cm")
		return cm;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'd':
	    if (token.Text()=="dyne")
		return dyne;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'e':
	    if (token.Text()=="erg")
		return erg;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'f':
	    if (token.Text()=="foot")
		return foot;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'g':
	    if ( u.size() == 1 )
		return g;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'i':
	    if (token.Text()=="inch")
		return inch;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'k':
	    if (token.Text()=="kg")
		return kg;
	    else if (token.Text()=="keV")
		return keV;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'l':
	    if (token.Text()=="lm")
		return lm;
	    else if (token.Text()=="lx")
		return lx;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'm':
	    if ( u.size() == 1 )
		return m;
	    else if (token.Text() == "mol")
		return mol;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'o':
	    if (token.Text()=="ohm")
		return ohm;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'p':
	    if (token.Text()=="pound")
		return pound;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 'r':
	    if (token.Text()=="rad")
		return rad;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	case 's':
	    if ( u.size() == 1 )
		return s;
	    else if (token.Text()=="sr")
		return sr;
	    else
		tokens.Report_Syntax_Error("expected a unit");

	default:
	    tokens.Report_Syntax_Error("expected a unit");
	}
    }
    else if (token.Type()==OTHER && token.Text()=="(")
    {
	Unit Result = Parse_Unit(tokens);
        token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!=")")
	    tokens.Report_Syntax_Error("missing ')'");
	return Result;
    }
    else
    {
	tokens.Report_Syntax_Error("expected a unit expression");
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

static Unit Parse_Unit_Term(Token_Stream &tokens)
{
    Unit const Result = Parse_Unit_Name(tokens);
    Token const token = tokens.Lookahead();
    if (token.Text()=="^")
    {
	tokens.Shift();
	double const exponent = Parse_Real(tokens);
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

Unit Parse_Unit(Token_Stream &tokens)
{
    if (!At_Unit_Term(tokens)) return dimensionless;

    Unit Result = Parse_Unit_Term(tokens);

    for (;;)
    {
	Token const token = tokens.Lookahead();
	if (token.Type()==OTHER)
	{
	    if (token.Text()=="-" && At_Unit_Term(tokens, 1))
	    {
		tokens.Shift();
		Result = Result * Parse_Unit_Term(tokens);
	    }
	    else if (token.Text()=="/")
	    {
		tokens.Shift();
		Result = Result / Parse_Unit_Term(tokens);
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
 * \brief Parse a dimensioned quantity.
 * 
 * \param tokens
 * Token stream from which to parse the quantity.
 * \param unit
 * Expected units for the quantity parsed, including conversion factor.
 * \param name
 * Name of the units expected for the quantity parsed, such as "length" or
 * "ergs/cm/sec/Hz". Used to generate diagnostic messages.
 *
 * \return The parsed value, converted to the desired unit system. 
 */

double Parse_Quantity(Token_Stream &tokens,
		      Unit const &target_unit,
		      char const *const name)
{
    double const value = Parse_Real(tokens);
    Unit const unit = Parse_Unit(tokens);
    if (!is_compatible(unit, target_unit))
    {
	ostringstream buffer;
	buffer << "expected quantity with dimensions of " << name;
	tokens.Report_Semantic_Error(buffer.str().c_str());
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

double Parse_Temperature(Token_Stream &tokens)
{
    double const T = Parse_Real(tokens);
    Unit const u = Parse_Unit(tokens);
    if (is_compatible(u, K))
    {
	double const Result = T * u.conv;
	if (Result<0.0)
	{
	    tokens.Report_Semantic_Error("temperature must be nonnegative");
	    return 0.0;
	}
	else
	{
	    return Result;
	}
    }
    else if (is_compatible(u, J))
    {
	double const Result = T * u.conv/rtt_units::boltzmannSI;
	if (Result<0.0)
	{
	    tokens.Report_Semantic_Error("energy must be nonnegative");
	    return 0.0;
	}
	else
	{
	    return Result;
	}
    }
    else
    {
	tokens.Report_Syntax_Error("expected quantity with units of "
				   "temperature");
	return 0.0;
    }
}

//---------------------------------------------------------------------------//
/*! 
 * Parses a STRING token and strips the delimiting quotation marks.
 *
 * \param tokens
 * Token_Stream from which to parse.
 * \return The stripped string.
 */

std::string Parse_Manifest_String(Token_Stream &tokens)
{
    Token const token = tokens.Shift();
    if (token.Type() != STRING)
    {
	tokens.Report_Syntax_Error("expected a string, but saw " + 
				   token.Text());
    }
    string Result = token.Text();
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

void Parse_Geometry(Token_Stream &tokens,
                    rtt_mesh_element::Geometry &parsed_geometry)
{
    if (parsed_geometry != rtt_mesh_element::END_GEOMETRY)
    {
        tokens.Report_Semantic_Error("geometry specified twice");
    }
    Token const token = tokens.Shift();
    if (token.Text() == "axisymmetric" ||
        token.Text() == "cylindrical")
    {
        parsed_geometry = rtt_mesh_element::AXISYMMETRIC;
    }
    else if (token.Text() == "cartesian" ||
             token.Text() == "xy" ||
             token.Text() == "slab")
    {
        parsed_geometry = rtt_mesh_element::CARTESIAN;
    }
    else if (token.Text() == "spherical")
	{
	    parsed_geometry = rtt_mesh_element::SPHERICAL;
	}
    else
    {
        tokens.Report_Syntax_Error(token,
                                   "expected a geometry option, but saw " + 
                                   token.Text());
    }
    Ensure(parsed_geometry == rtt_mesh_element::AXISYMMETRIC ||
           parsed_geometry == rtt_mesh_element::CARTESIAN    ||
           parsed_geometry == rtt_mesh_element::SPHERICAL);
    return;
}

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of utilities.cc
//---------------------------------------------------------------------------//
