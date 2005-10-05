//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstConsole_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed May 19 11:26:15 MDT 2004
 * \brief  Unit tests for Console_Token_Stream class.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>

#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include "parser_test.hh"
#include "../Release.hh"
#include "../Console_Token_Stream.hh"

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstConsole_Token_Stream()
{
    {
	set<char> ws;
	ws.insert(' ');
	Console_Token_Stream tokens(ws);
	if (tokens.Whitespace()!=ws)
	{
	    FAILMSG("whitespace characters are NOT correctly set");
	}
	else
	{
	    PASSMSG("whitespace characters are correctly set");
	}	
    }

    {
	Console_Token_Stream tokens;
	if (tokens.Whitespace()!=Text_Token_Stream::default_whitespace)
	{
	    FAILMSG("whitespace characters are NOT correct defaults");
	}
	else
	{
	    PASSMSG("whitespace characters are correct defaults");
	}

	Token token = tokens.Lookahead(3);
	if (token.Type()!=KEYWORD || token.Text()!="COLOR") 
	{
	    FAILMSG("Lookahead(3) does NOT have correct value");
	}
	else
	{
	    PASSMSG("Lookahead(3) has correct value");
	}

	tokens.Report_Semantic_Error(token, "dummy error");
	if (tokens.Error_Count()!=1)
	{
	    FAILMSG("Dummy error NOT counted properly");
	}
	else
	{
	    PASSMSG("Dummy error counted properly");
	}

	tokens.Report_Semantic_Error("dummy error");
	if (tokens.Error_Count()!=2)
	{
	    FAILMSG("Dummy error NOT counted properly");
	}
	else
	{
	    PASSMSG("Dummy error counted properly");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="BLUE")
	{
	    FAILMSG("First shift does NOT have correct value");
	}
	else
	{
	    PASSMSG("First shift has correct value");
	}

	token = tokens.Lookahead();
	if (token.Type()!=KEYWORD || token.Text()!="GENERATE ERROR")
	{
	    FAILMSG("Lookahed after first shift does NOT have correct value");
	}
	else
	{
	    PASSMSG("Lookahead after first shift has correct value");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="GENERATE ERROR")
	{
	    FAILMSG("Second shift does NOT have correct value");
	}
	else
	{
	    PASSMSG("Second shift has correct value");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || 
	    token.Text()!="GENERATE ANOTHER ERROR")
	{
	    FAILMSG("Third shift does NOT have correct value");
	}
	else
	{
	    PASSMSG("Third shift has correct value");
	}

        token = Token('$', "test_parser");
	tokens.Pushback(token);

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="$")
	{
	    FAILMSG("Shift after pushback does NOT have correct value");
	}
	else
	{
	    PASSMSG("Shift after pushback has correct value");
	}

	try 
	{
	    tokens.Report_Syntax_Error(token, "dummy syntax error");  
	    FAILMSG("Syntax error NOT correctly thrown");
	}
	catch (const Syntax_Error &msg)
	{
	    PASSMSG("Syntax error correctly thrown and caught");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="COLOR") ITFAILS;
	
	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="BLACK") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=END) ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="-1.563e+3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.563e+3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!=".563e+3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!=".") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="-") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.563") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.e+3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.e3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1e+3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1e3") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="19090") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="01723") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="0x1111a") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="0") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="8123") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=STRING || token.Text()!="\"manifest string\"")
	    ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=STRING || 
	    token.Text()!="\"manifest \\\"string\\\"\"")
	    ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="@") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="1") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="e") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=INTEGER || token.Text()!="0") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="x") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=EXIT) ITFAILS;
	token = tokens.Shift();
	if (token.Type()!=EXIT) ITFAILS;

	tokens.Rewind();
	token = tokens.Lookahead();
	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="BLUE") ITFAILS;
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
	if (rtt_c4::nodes() == 1) tstConsole_Token_Stream();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstConsole_Token_Stream, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_parser_test::passed) 
    {
        cout << "**** tstConsole_Token_Stream Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    rtt_c4::global_barrier();
    cout << "Done testing tstConsole_Token_Stream." << endl;
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstConsole_Token_Stream.cc
//---------------------------------------------------------------------------//
