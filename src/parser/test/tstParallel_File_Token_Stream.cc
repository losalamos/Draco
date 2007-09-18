//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstParallel_File_Token_Stream.cc
 * \author Kent Budge
 * \date   Fri Apr  4 09:34:28 2003
 * \brief  Unit tests for class Parallel_File_Token_Stream
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"

#include "parser_test.hh"
#include "../Release.hh"
#include "../Parallel_File_Token_Stream.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstParallel_File_Token_Stream()
{
    using namespace rtt_parser;

    {
	Parallel_File_Token_Stream tokens("scanner_test.inp");
	if (tokens.Whitespace()!=Text_Token_Stream::default_whitespace)
	{
	    FAILMSG("Whitespace not set correctly");
	}
	else
	{
	    PASSMSG("Whitespace set correctly.");
	}

	Token token = tokens.Lookahead(4);
	if (token.Type()!=KEYWORD || token.Text()!="BLACK")
	{
	    FAILMSG("Keyword not read correctly");
	}
	else
	{
	    PASSMSG("Keyword read correctly.");
	}

	tokens.Report_Semantic_Error(token, "dummy error");
	if (tokens.Error_Count()!=1)
	{
	    FAILMSG("Semantic error not handled correctly.");
	}
	else
	{
	    PASSMSG("Semantic error handled correctly.");
	}

	tokens.Report_Semantic_Error("dummy error");
	if (tokens.Error_Count()!=2)
	{
	    FAILMSG("Second semantic error not handled correctly.");
	}
	else
	{
	    PASSMSG("Second semantic error handled correctly.");
	}
    }

    {
	set<char> ws;
	ws.insert(':');
	Parallel_File_Token_Stream tokens("scanner_test.inp", ws);
	if (tokens.Whitespace()!=ws)
	{
	    FAILMSG("Whitespace not set correctly");
	}
	else
	{
	    PASSMSG("Whitespace set correctly.");
	}

	Token token = tokens.Lookahead(4);
	if (token.Type()!=OTHER || token.Text()!="=")
	{
	    FAILMSG("'=' token not read correctly");
	}
	else
	{
	    PASSMSG("'=' token read correctly.");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="BLUE")
	{
	    FAILMSG("Keyword BLUE not read correctly");
	}
	else
	{
	    PASSMSG("Keyword BLUE read correctly.");
	}
	token = tokens.Lookahead();
	if (token.Type()!=KEYWORD || token.Text()!="GENERATE ERROR")
	{
	    FAILMSG("Keyword GENERATE ERROR not read correctly");
	}
	else
	{
	    PASSMSG("Keyword GENERATE ERROR read correctly.");
	}

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || token.Text()!="GENERATE ERROR") ITFAILS;

	token = tokens.Shift();
	if (token.Type()!=KEYWORD || 
	    token.Text()!="GENERATE ANOTHER ERROR") ITFAILS;

        token = Token('$', "test_parser");
	tokens.Pushback(token);

	token = tokens.Shift();
	if (token.Type()!=OTHER || token.Text()!="$") ITFAILS;

	try 
	{
	    tokens.Report_Syntax_Error(token, "dummy syntax error");  
	    {
		ostringstream msg;
		msg << "Parallel_File_Token_Stream did not throw an exception when\n"
		    << "\ta syntax error was reported by Token_Stream." << endl;
		FAILMSG( msg.str() );
	    }
	}
	catch (const Syntax_Error &)
	{
	    {
		ostringstream msg;
		msg << "Parallel_File_Token_Stream threw an expected exception when\n"
		    << "\ta syntax error was reported by Token_Stream." << endl;
		PASSMSG( msg.str() );
	    }
	}
	if (tokens.Error_Count()!=1) ITFAILS;

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

        // Check invariance even when --with-dbc=0.
        if( ! tokens.check_class_invariants() ) ITFAILS;
    }

    {
	try
	{
	    Parallel_File_Token_Stream tokens("no such file");
	    // The preceeding file does not exist.
	    ostringstream errmsg;
	    errmsg << "Parallel_File_Token_Stream did not throw an expected exception.\n"
		   << "\tThe constructor should throw an exception if the requested\n"
		   << "\tfile can not be opened." << endl;
	    FAILMSG( errmsg.str() );
	    // Token token = tokens.Shift();
	    // if (token.Type()!=ERROR) ITFAILS;
	}
	catch( std::invalid_argument const &a )
	{
	    std::ostringstream errmsg;
	    errmsg << "Parallel_File_Token_Stream threw an expected exception.\n"
		   << "\tThe constructor should throw an exception if the requested\n"
		   << "\tfile can not be opened." << endl;
	    PASSMSG( errmsg.str() );
	}
	catch ( ... )
	{
	    ostringstream errmsg;
	    errmsg << "Parallel_File_Token_Stream threw an unknown exception "
		   << "during contruction." << endl;
	    FAILMSG( errmsg.str() );
	}
    }

    {
	Parallel_File_Token_Stream tokens("scanner_recovery.inp");
	bool exception = false;
	try
	{
	    tokens.Shift();
	}
	catch (const Syntax_Error &msg)
	{
	    cout << msg.what() << endl;
	    exception = true;
	}
	if (!exception) ITFAILS;

	exception = false;
	try
	{
	    tokens.Shift();
	}
	catch (const Syntax_Error &msg)
	{
	    cout << msg.what() << endl;
	    exception = true;
	}
	if (!exception) ITFAILS;
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
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_parser::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	tstParallel_File_Token_Stream();
    }
    catch (std::exception &ass)
    {
	cout << "While testing tstParallel_File_Token_Stream, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }
    catch (...)
    {
	cout << "While testing tstParallel_File_Token_Stream, "
             "an unrecognized exception was thrown."
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_parser_test::passed) 
	{
	    cout << "**** tstParallel_File_Token_Stream Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing tstParallel_File_Token_Stream on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParallel_File_Token_Stream.cc
//---------------------------------------------------------------------------//
