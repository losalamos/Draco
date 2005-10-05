//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstParse_Table.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for the Parse_Table class.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original revision
 * 1) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "parser_test.hh"
#include "../Release.hh"
#include "../Parse_Table.hh"
#include "../File_Token_Stream.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

using namespace std;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

static const char *color[2] = {"BLACK", "BLUE"};
bool color_set[2];

static void Parse_Color(Token_Stream &, int i)
{
  cout << "You have requested " << color[i] << endl;
  color_set[i] = true;
}

static void Parse_Any_Color(Token_Stream &tokens, int)
{
  Token token = tokens.Shift();
  for (unsigned i=0; i<sizeof(color)/sizeof(const char*); i++)
    if (!strcmp(token.Text().c_str(), color[i])){
      cout << "You have requested " << color[i] << endl;
      color_set[i] = true;
      return;
    }

  tokens.Report_Syntax_Error(token, "expected a color");
}     

const Keyword raw_table[] = {
  {"BLUE", Parse_Color, 1, "main"},
  {"BLACK", Parse_Color, 0, "main"},
  {"COLOR", Parse_Any_Color, 0, "main"},
};
const size_t raw_table_size = sizeof(raw_table)/sizeof(Keyword);

const Keyword raw_table_2[] = {
  {"BLUE", Parse_Color, 1, "main"},
  {"BLACK", Parse_Color, 0, "main"},
};
const size_t raw_table_2_size = sizeof(raw_table_2)/sizeof(Keyword);

void tstKeyword()
{
    Keyword black = {"BLACK", Parse_Color, 0, "main"};

    if (black == raw_table_2[1])
    {
	PASSMSG("keyword equality operator okay");
    }
    else
    {
	FAILMSG("keyword equality operator NOT okay");
    }
}

void tstParse_Table()
{
    Parse_Table table;
    
    table.reserve(raw_table_size);
    table.Add(raw_table, raw_table_size);

    if (table.size()!=raw_table_size) ITFAILS;
    
    File_Token_Stream token_stream("parser_test.inp");
    
    table.Parse(token_stream);
    
    if (!color_set[1]) ITFAILS;
    
    if (token_stream.Error_Count()!=5) ITFAILS;


    token_stream.Rewind();

    table.Set_Flags(Parse_Table::CASE_INSENSITIVE);

    color_set[0] = color_set[1] = 0;
    table.Parse(token_stream);

    if (!color_set[1]) ITFAILS;
    if (token_stream.Error_Count()!=4) ITFAILS;


    token_stream.Rewind();

    table.Set_Flags(Parse_Table::CASE_INSENSITIVE |
		    Parse_Table::PARTIAL_IDENTIFIER_MATCH);

    color_set[0] = color_set[1] = 0;
    table.Parse(token_stream);

    if (!color_set[1]) ITFAILS;
    if (token_stream.Error_Count()!=3) ITFAILS;

    Parse_Table table_2(raw_table, raw_table_size);

    if (table_2.size()!=raw_table_size) ITFAILS;
    
    token_stream.Rewind();
    
    table_2.Parse(token_stream);
    
    if (!color_set[1]) ITFAILS;
    
    if (token_stream.Error_Count()!=5) ITFAILS;

    Keyword test_key = {"THIS SHOULD WORK", Parse_Color, 0, 0};
    if (!Is_Well_Formed_Keyword(test_key)) ITFAILS;

    Keyword benign_ambiguous_table[] = 
	{
	    {"KEY", Parse_Color, 0, 0},
	    {"KEY", Parse_Color, 0, 0}
	};
    table_2.Add(benign_ambiguous_table, 2);
    token_stream.Rewind();
    table_2.Parse(token_stream);

    Keyword malign_ambiguous_table[] = 
	{
	    {"KEY", Parse_Color, 1, 0}
	};
    try 
    {
	table_2.Add(malign_ambiguous_table, 1);
	token_stream.Rewind();
	table_2.Parse(token_stream);
	FAILMSG("did NOT catch ambiguous keyword");
    }
    catch (const rtt_dsxx::assertion &msg)
    {
	cout << msg.what() << endl;
	PASSMSG("successfully detected ambiguous keyword");
    }

    File_Token_Stream recover_stream("recovery.inp");
    table.Parse(recover_stream);
    if (recover_stream.Error_Count() != 2) ITFAILS;

    Parse_Table table_3;
    Keyword case_ambiguous_table[] = 
	{
	    {"key", Parse_Color, 0, 0},
	    {"Key", Parse_Color, 1, 0}
	};    
    try 
    {
	table_3.Add(case_ambiguous_table, 2);
	table_3.Parse(token_stream);
	table_3.Set_Flags(Parse_Table::CASE_INSENSITIVE);
	token_stream.Rewind();
	table_3.Parse(token_stream);
	FAILMSG("did NOT catch case-dependent ambiguous keyword");
    }
    catch (const rtt_dsxx::assertion &msg)
    {
	cout << msg.what() << endl;
	PASSMSG("successfully detected case-dependent ambiguous keyword");
    }

    {
	// Additional test mandated by bug discovery.

	Parse_Table table;
	
	table.reserve(raw_table_2_size);
	table.Add(raw_table_2, raw_table_2_size);
	
	if (table.size()!=raw_table_2_size) ITFAILS;
	
	File_Token_Stream token_stream("parser_test.inp");
	
	table.Parse(token_stream);

	if (token_stream.Error_Count() != 5) ITFAILS;
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
	if ( rtt_c4::nodes() == 1 )
	{
	    tstKeyword();
	    tstParse_Table();
	}
    }
    catch (std::exception &ass)
    {
	cout << "While testing tstParse_Table, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_parser_test::passed) 
    {
        cout << "**** tstParse_Table Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;
    rtt_c4::global_barrier();
    cout << "Done testing tstParse_Table." << endl;
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParse_Table.cc
//---------------------------------------------------------------------------//
