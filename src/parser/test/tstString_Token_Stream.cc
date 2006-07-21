//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstString_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for String_Token_Stream class.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>

#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include "parser_test.hh"
#include "../Release.hh"
#include "../String_Token_Stream.hh"

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstString_Token_Stream()
{
    ifstream infile("scanner_test.inp");
    string contents;
    while (true)
    {
	char const c = infile.get();
	if (infile.eof() || infile.fail()) break;
	contents += c;
    }

    {
	String_Token_Stream tokens(contents);
	if (tokens.Whitespace()!=Text_Token_Stream::default_whitespace)
	{
	    FAILMSG("whitespace characters are NOT correct defaults");
	}
	else
	{
	    PASSMSG("whitespace characters are correct defaults");
	}

	Token token = tokens.Lookahead(4);
	if (token.Type()!=KEYWORD || token.Text()!="BLACK") 
	{
	    FAILMSG("Lookahead(4) does NOT have correct value");
	}
	else
	{
	    PASSMSG("Lookahead(4) has correct value");
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

        if( ! tokens.check_class_invariants() ) ITFAILS;
    }

    {
	set<char> ws;
	ws.insert(':');
	String_Token_Stream tokens(contents, ws);
	if (tokens.Whitespace()!=ws)
	{
	    FAILMSG("whitespace characters are NOT correctly specified");
	}
	else
	{
	    PASSMSG("whitespace characters are correctly specified");
	}

	Token token = tokens.Lookahead(4);
	if (token.Type()!=OTHER || token.Text()!="=")
	{
	    FAILMSG("Lookahead(4) does NOT have correct value");
	}
	else
	{
	    PASSMSG("Lookahead(4) has correct value");
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
	if (tokens.Error_Count()!=1)
	{
	    FAILMSG("Syntax error NOT correctly counted");
	}
	else
	{
	    PASSMSG("Syntax error correctly counted");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="COLOR") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="=") ITFAILS;
	
	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="BLACK") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=END) ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="-") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=REAL || token.Text()!="1.563e+3") ITFAILS;

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

//---------------------------------------------------------------------------//

    {
	ifstream infile("scanner_recovery.inp");
	string contents;
	while (true)
	{
	    char const c = infile.get();
	    if (infile.eof() || infile.fail()) break;
	    contents += c;
	}
	String_Token_Stream tokens(contents);
	try
	{
	    tokens.Shift();
	    ostringstream msg;
	    msg << "Token_Stream did not throw an exception when\n"
		<< "\tunbalanced quotes were read from the input\n"
		<< "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
	    FAILMSG( msg.str() );
	}
	catch ( const Syntax_Error &msg )
	{
	    // cout << msg.what() << endl;
	    // exception = true;
	    string errmsg = msg.what();
	    string expected( "syntax error" );
	    if( errmsg == expected )
	    {
		ostringstream msg;
		msg << "Caught expected exception from Token_Stream.\n"
		    << "\tunbalanced quotes were read from the input\n"
		    << "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
		PASSMSG( msg.str() );
	    }
	    else ITFAILS;
	}

	try
	{
	    tokens.Shift();
	    ostringstream msg;
	    msg << "Token_Stream did not throw an exception when\n"
		<< "\tunbalanced quotes were read from the input\n"
		<< "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
	    FAILMSG( msg.str() );
	}
	catch  (const Syntax_Error &msg )
	{
	    //cout << msg.what() << endl;
	    // exception = true;
	    string errmsg = msg.what();
	    string expected( "syntax error" );
	    if( errmsg == expected )
	    {
		ostringstream msg;
		msg << "Caught expected exception from Token_Stream.\n"
		    << "\tunbalanced quotes were read from the input\n"
		    << "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
		PASSMSG( msg.str() );
	    }
	    else ITFAILS;
	}
    }

    return;
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
	if (rtt_c4::nodes() == 1) tstString_Token_Stream();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstString_Token_Stream, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_parser_test::passed) 
    {
        cout << "**** tstString_Token_Stream Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    rtt_c4::global_barrier();
    cout << "Done testing tstString_Token_Stream." << endl;
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstString_Token_Stream.cc
//---------------------------------------------------------------------------//
