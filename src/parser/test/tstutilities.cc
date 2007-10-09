//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstutilities.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <limits>
#include "ds++/Soft_Equivalence.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser_test.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../File_Token_Stream.hh"
#include "../String_Token_Stream.hh"
#include "../utilities.hh"
#include "../Unit.hh"

using namespace std;

using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstutilities(UnitTest &ut)
{
    File_Token_Stream tokens("utilities.inp");

    // Try to read a real number.

    double d = parse_real(tokens);
    if (tokens.error_count() != 0 || d != 5.)
    {
	FAILMSG("real NOT successfully parsed");
    }
    else
    {
	PASSMSG("real successfully parsed");
    }

    // Try to read an integer.

    int i = parse_integer(tokens);
    if (tokens.error_count() != 0 || i != 1)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read a negative integer.

    i = parse_integer(tokens);
    if (tokens.error_count() != 0 || i != -3)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read an unsigned integer.

    i = parse_unsigned_integer(tokens);
    if (tokens.error_count() != 0 || i != 4)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read a positive integer.

    i = parse_positive_integer(tokens);
    if (tokens.error_count() != 0 || i != 1198)
    {
	FAILMSG("positive integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("positive integer successfully parsed");
    }

    // Try to read an integer as a real.

    d = parse_real(tokens);
    if (tokens.error_count() != 0 || d != 2.)
    {
	FAILMSG("integer NOT successfully parsed as real");
    }
    else
    {
	PASSMSG("integer successfully parsed as real");
    }

    // Try to read some vectors.

    double v[3];
    parse_vector(tokens, v);
    Token token = tokens.shift();
    if (v[0] == 3. && v[1] == 0.0 && v[2] == 0.0 && 
	token.type() == KEYWORD && token.text() == "stop")
    {
	PASSMSG("1-D vector successfully parsed");
    }
    else
    {
	FAILMSG("1-D vector NOT successfully parsed");
    }

    parse_vector(tokens, v);
    token = tokens.shift();
    if (v[0] == 1. && v[1] == 2.0 && v[2] == 0.0 && 
	token.type() == KEYWORD && token.text() == "stop")
    {
	PASSMSG("2-D vector successfully parsed");
    }
    else
    {
	FAILMSG("2-D vector NOT successfully parsed");
    }

    parse_vector(tokens, v);
    if (v[0] == 4. && v[1] == 3.0 && v[2] == 2.0
	&& tokens.shift().text()=="stop")
    {
	PASSMSG("3-D vector successfully parsed");
    }
    else
    {
	FAILMSG("3-D vector NOT successfully parsed");
    }
    unsigned w[3];
    parse_unsigned_vector(tokens, w, 3);
    token = tokens.shift();
    if (w[0] == 3 && w[1] == 2 && w[2] == 1 &&
        token.type() == KEYWORD && token.text() == "stop")
    {
        PASSMSG("Vector of unsigned successfully parsed");
    }
    else
    {
        FAILMSG("Vector of unsigned NOT successfully parsed");
    }
    
    // Try to read some unit expressions

    Unit one = {0,0,0,0,0,0,0,0,0, 1};

    Unit left = parse_unit(tokens);
    if (left!=J) ITFAILS;

    left = parse_unit(tokens);
    Unit right = parse_unit(tokens);
    if (left!=right || left!=C) ITFAILS;

    left = parse_unit(tokens);
    cout << left << endl;
    if (left!=1/s) ITFAILS;

    left = parse_unit(tokens);
    right = parse_unit(tokens);
    if (left!=right || left!=N) ITFAILS;

    left = parse_unit(tokens);
    if (left!=J/K) ITFAILS;

    left = parse_unit(tokens);
    if (left!=J/mol) ITFAILS;

    left = parse_unit(tokens);
    right = parse_unit(tokens);
    if (left!=right || left!=lm) ITFAILS;

    left = parse_unit(tokens);
    if (left!=rad/s) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one) ITFAILS;

    left = parse_unit(tokens);
    if (left!=one)
    {
	FAILMSG("dyne definition did NOT check out");
    }
    else
    {
	PASSMSG("dyne definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=one)
    {
	FAILMSG("erg definition did NOT check out");
    }
    else
    {
	PASSMSG("erg definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, cm) || !soft_equiv(left.conv, 0.0254))
    {
	FAILMSG("inch definition did NOT check out");
    }
    else
    {
	PASSMSG("inch definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 12.0))
    {
	FAILMSG("foot definition did NOT check out");
    }
    else
    {
	PASSMSG("foot definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 4.448221615))
    {
	FAILMSG("pound definition did NOT check out");
    }
    else
    {
	PASSMSG("pound definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one))
    {
	FAILMSG("keV definition did NOT check out");
    }
    else
    {
	PASSMSG("keV definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=J) ITFAILS;

    left = parse_unit(tokens);
    if (left!=J) ITFAILS;

    left = parse_unit(tokens);
    if (left!=K)
    {
	FAILMSG("K definition did NOT check out");
    }
    else
    {
	PASSMSG("K definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=sr)
    {
	FAILMSG("sr definition did NOT check out");
    }
    else
    {
	PASSMSG("sr definition checks out");
    }

    // Now see if we catch a bogus unit expression.
    try
    {
	left = parse_unit(tokens);
	FAILMSG("did NOT successfully catch bogus unit");
    }
    catch (const Syntax_Error &)
    {
	PASSMSG("successfully caught bogus unit");
    }



    // Try to read some dimensioned quantities.
    double length = parse_quantity(tokens, rtt_parser::m, "length");
    if (fabs(length-3.0)<=std::numeric_limits<double>::epsilon())
    {
	PASSMSG("length successfully parsed");
    }
    else
    {
	FAILMSG("length NOT successfully parsed");
    }

    double energy = parse_quantity(tokens, rtt_parser::J, "energy");
    if (fabs(energy-2.3e-7)<=std::numeric_limits<double>::epsilon())
    {
	PASSMSG("cgs energy successfully parsed");
    }
    else
    {
	FAILMSG("cgs energy NOT successfully parsed");
    }

    unsigned old_error_count = tokens.error_count();
    length = parse_quantity(tokens, rtt_parser::m, "length");
    if (tokens.error_count()==old_error_count)
    {
	FAILMSG("bad length NOT successfully detected");
    }
    else
    {
	PASSMSG("bad length successfully detected");
    }

    old_error_count = tokens.error_count();
    double T = parse_temperature(tokens);
    if (tokens.error_count()!=old_error_count)
    {
	FAILMSG("temperature NOT successfully parsed");
    }
    else
    {
	PASSMSG("temperature successfully parsed");
    }
    T = parse_temperature(tokens);
    if (tokens.error_count()!=old_error_count)
    {
	FAILMSG("temperature NOT successfully parsed");
    }
    else
    {
	PASSMSG("temperature successfully parsed");
    }

    // Try reading sequence of quantities with signs
    T = parse_quantity(tokens, J, "energy");
    if (tokens.error_count()!=old_error_count)
    {
	FAILMSG("second negative quantity NOT successfully parsed");
    }
    else
    {
	PASSMSG("second negative quantity successfully parsed");
    }

    // Try reading a manifest string.
    string parsed_string = parse_manifest_string(tokens);
    if (parsed_string!="manifest string")
    {
	FAILMSG("manifest string NOT successfully parsed");
    }
    else
    {
	PASSMSG("manifest string successfully parsed");
    }

    // Try reading a geometry.
    rtt_mesh_element::Geometry geometry = rtt_mesh_element::END_GEOMETRY;
    parse_geometry(tokens, geometry);
    if (geometry != rtt_mesh_element::AXISYMMETRIC)
    {
	FAILMSG("geometry NOT successfully parsed");
    }
    else
    {
	PASSMSG("geometry successfully parsed");
    }

    String_Token_Stream string("4.5");
    if (soft_equiv(parse_positive_real(string), 4.5))
    {
        ut.passes("read positive real correctly");
    }
    else
    {
        ut.failure("did NOT read positive real correctly");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);

    try
    {
	// >>> UNIT TESTS
        tstutilities(ut);
    }
    catch (std::exception &ass)
    {
	cout << "While testing tstutilities, " << ass.what()
	     << endl;
	ut.numFails++;
    }
    catch (...)
    {
        cout << "ERROR: While testing tstBudge_Opacity_Model, "
             << "An unknown exception was thrown." << endl;
        ut.numFails++;
    }

    if (ut.numFails==0)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}   

//---------------------------------------------------------------------------//
//                        end of tstutilities.cc
//---------------------------------------------------------------------------//
