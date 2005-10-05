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
#include "parser_test.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../File_Token_Stream.hh"
#include "../utilities.hh"
#include "../Unit.hh"

using namespace std;

using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstutilities()
{
    File_Token_Stream tokens("utilities.inp");

    // Try to read a real number.

    double d = Parse_Real(tokens);
    if (tokens.Error_Count() != 0 || d != 5.)
    {
	FAILMSG("real NOT successfully parsed");
    }
    else
    {
	PASSMSG("real successfully parsed");
    }

    // Try to read an integer.

    int i = Parse_Integer(tokens);
    if (tokens.Error_Count() != 0 || i != 1)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read a negative integer.

    i = Parse_Integer(tokens);
    if (tokens.Error_Count() != 0 || i != -3)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read an unsigned integer.

    i = Parse_Unsigned_Integer(tokens);
    if (tokens.Error_Count() != 0 || i != 4)
    {
	FAILMSG("integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("integer successfully parsed");
    }

    // Try to read a positive integer.

    i = Parse_Positive_Integer(tokens);
    if (tokens.Error_Count() != 0 || i != 1198)
    {
	FAILMSG("positive integer NOT successfully parsed");
    }
    else
    {
	PASSMSG("positive integer successfully parsed");
    }

    // Try to read an integer as a real.

    d = Parse_Real(tokens);
    if (tokens.Error_Count() != 0 || d != 2.)
    {
	FAILMSG("integer NOT successfully parsed as real");
    }
    else
    {
	PASSMSG("integer successfully parsed as real");
    }

    // Try to read some vectors.

    double v[3];
    Parse_Vector(tokens, v);
    Token token = tokens.Shift();
    if (v[0] == 3. && v[1] == 0.0 && v[2] == 0.0 && 
	token.Type() == KEYWORD && token.Text() == "stop")
    {
	PASSMSG("1-D vector successfully parsed");
    }
    else
    {
	FAILMSG("1-D vector NOT successfully parsed");
    }

    Parse_Vector(tokens, v);
    token = tokens.Shift();
    if (v[0] == 1. && v[1] == 2.0 && v[2] == 0.0 && 
	token.Type() == KEYWORD && token.Text() == "stop")
    {
	PASSMSG("2-D vector successfully parsed");
    }
    else
    {
	FAILMSG("2-D vector NOT successfully parsed");
    }

    Parse_Vector(tokens, v);
    if (v[0] == 4. && v[1] == 3.0 && v[2] == 2.0
	&& tokens.Shift().Text()=="stop")
    {
	PASSMSG("3-D vector successfully parsed");
    }
    else
    {
	FAILMSG("3-D vector NOT successfully parsed");
    }

    // Try to read some unit expressions

    Unit one = {0,0,0,0,0,0,0,0,0, 1};

    Unit left = Parse_Unit(tokens);
    if (left!=J) ITFAILS;

    left = Parse_Unit(tokens);
    Unit right = Parse_Unit(tokens);
    if (left!=right || left!=C) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=1/s) ITFAILS;

    left = Parse_Unit(tokens);
    right = Parse_Unit(tokens);
    if (left!=right || left!=N) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=J/K) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=J/mol) ITFAILS;

    left = Parse_Unit(tokens);
    right = Parse_Unit(tokens);
    if (left!=right || left!=lm) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=rad/s) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=one)
    {
	FAILMSG("dyne definition did NOT check out");
    }
    else
    {
	PASSMSG("dyne definition checks out");
    }

    left = Parse_Unit(tokens);
    if (left!=one)
    {
	FAILMSG("erg definition did NOT check out");
    }
    else
    {
	PASSMSG("erg definition checks out");
    }

    left = Parse_Unit(tokens);
    if (!is_compatible(left, cm) || !soft_equiv(left.conv, 0.0254))
    {
	FAILMSG("inch definition did NOT check out");
    }
    else
    {
	PASSMSG("inch definition checks out");
    }

    left = Parse_Unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 12.0))
    {
	FAILMSG("foot definition did NOT check out");
    }
    else
    {
	PASSMSG("foot definition checks out");
    }

    left = Parse_Unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 4.448221615))
    {
	FAILMSG("pound definition did NOT check out");
    }
    else
    {
	PASSMSG("pound definition checks out");
    }

    left = Parse_Unit(tokens);
    if (!is_compatible(left, one))
    {
	FAILMSG("keV definition did NOT check out");
    }
    else
    {
	PASSMSG("keV definition checks out");
    }

    left = Parse_Unit(tokens);
    if (left!=J) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=J) ITFAILS;

    left = Parse_Unit(tokens);
    if (left!=K)
    {
	FAILMSG("K definition did NOT check out");
    }
    else
    {
	PASSMSG("K definition checks out");
    }

    left = Parse_Unit(tokens);
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
	left = Parse_Unit(tokens);
	FAILMSG("did NOT successfully catch bogus unit");
    }
    catch (const Syntax_Error &)
    {
	PASSMSG("successfully caught bogus unit");
    }



    // Try to read some dimensioned quantities.
    double length = Parse_Quantity(tokens, rtt_parser::m, "length");
    if (fabs(length-3.0)<=std::numeric_limits<double>::epsilon())
    {
	PASSMSG("length successfully parsed");
    }
    else
    {
	FAILMSG("length NOT successfully parsed");
    }

    double energy = Parse_Quantity(tokens, rtt_parser::J, "energy");
    if (fabs(energy-2.3e-7)<=std::numeric_limits<double>::epsilon())
    {
	PASSMSG("cgs energy successfully parsed");
    }
    else
    {
	FAILMSG("cgs energy NOT successfully parsed");
    }

    unsigned old_error_count = tokens.Error_Count();
    length = Parse_Quantity(tokens, rtt_parser::m, "length");
    if (tokens.Error_Count()==old_error_count)
    {
	FAILMSG("bad length NOT successfully detected");
    }
    else
    {
	PASSMSG("bad length successfully detected");
    }

    old_error_count = tokens.Error_Count();
    double T = Parse_Temperature(tokens);
    if (tokens.Error_Count()!=old_error_count)
    {
	FAILMSG("temperature NOT successfully parsed");
    }
    else
    {
	PASSMSG("temperature successfully parsed");
    }
    T = Parse_Temperature(tokens);
    if (tokens.Error_Count()!=old_error_count)
    {
	FAILMSG("temperature NOT successfully parsed");
    }
    else
    {
	PASSMSG("temperature successfully parsed");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_parser::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	if (rtt_c4::nodes() == 1) tstutilities();
    }
    catch (std::exception &ass)
    {
	cout << "While testing tstutilities, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_parser_test::passed) 
    {
        cout << "**** tstutilities Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    rtt_c4::global_barrier();
    cout << "Done testing tstutilities." << endl;
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstutilities.cc
//---------------------------------------------------------------------------//
