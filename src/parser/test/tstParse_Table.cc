//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstParse_Table.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for the Parse_Table class.
 * \note   Copyright � 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original revision
 * 1) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <string.h>
#include "parser_test.hh"
#include "../Release.hh"
#include "../Parse_Table.hh"
#include "../File_Token_Stream.hh"
#include "../String_Token_Stream.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

static const char *color[3] = {"BLACK", "BLUE", "BLUE GREEN"};
bool color_set[3];

static void Parse_Color(Token_Stream &, int i)
{
  cout << "You have requested " << color[i] << endl;
  color_set[i] = true;
}

static void Parse_Any_Color(Token_Stream &tokens, int)
{
  Token token = tokens.shift();
  for (unsigned i=0; i<sizeof(color)/sizeof(const char*); i++)
    if (!strcmp(token.text().c_str(), color[i])){
      cout << "You have requested " << color[i] << endl;
      color_set[i] = true;
      return;
    }

  tokens.report_syntax_error(token, "expected a color");
}     

const Keyword raw_table[] = {
  {"BLUE", Parse_Color, 1, "main"},
  {"BLACK", Parse_Color, 0, "main"},
  {"BLUE GREEN", Parse_Color, 2, "main"},
  {"BLUISH GREEN", Parse_Color, 2, "main"},
  {"lower blue", Parse_Color, 2, "main"},
  {"COLOR", Parse_Any_Color, 0, "main"},
};
const size_t raw_table_size = sizeof(raw_table)/sizeof(Keyword);

const Keyword raw_table_2[] = {
  {"BLUE", Parse_Color, 1, "main"},
  {"BLACK", Parse_Color, 0, "main"},
};
const size_t raw_table_2_size = sizeof(raw_table_2)/sizeof(Keyword);

class Error_Token_Stream : public Token_Stream
{
  public:

    void rewind(){}

  protected:

    void report(Token const &, string const &err)
    {
        cout << "error reported to Error_Token_Stream" << endl;
    }

    void report(string const &err)
    {
        cout << "error reported to Error_Token_Stream" << endl;
    }

    Token fill_()
    {
        return Token(ERROR, "error");
    }
};

class Colon_Token_Stream : public Token_Stream
{
  public:

    Colon_Token_Stream() : count_(0) {}

    void rewind(){}

  protected:

    void report(Token const &, string const &err)
    {
        cout << "error reported to Colon_Token_Stream" << endl;
    }

    void report(string const &err)
    {
        cout << "error reported to Colon_Token_Stream" << endl;
    }

    Token fill_()
    {
        switch (count_++)
        {
            case 0:
                return Token(';', "");
            case 1:
                return Token(END, "end");
            case 2:
                return Token(EXIT, "");
            default:
                Insist(false, "bad case");
                return Token(ERROR, "");
        }
    }

  private:

    unsigned count_;
};

//---------------------------------------------------------------------------//
void tstKeyword(UnitTest &ut)
{
    Keyword black = {"BLACK", Parse_Color, 0, "main"};

    if (black == raw_table_2[1])
    {
	ut.passes("keyword equality operator okay");
    }
    else
    {
	ut.failure("keyword equality operator NOT okay");
    }

    {
        Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
        Keyword key2 = {"RED", Parse_Color, 1, "main"};
        if (key1==key2) ut.failure("comparison of dislike keywords FAILED");
    }
    {
        Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
        Keyword key2 = {"BLUE", Parse_Any_Color, 1, "main"};
        if (key1==key2) ut.failure("comparison of dislike functions FAILED");
    }
    {
        Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
        Keyword key2 = {"BLUE", Parse_Color, 12, "main"};
        if (key1==key2) ut.failure("comparison of dislike index FAILED");
    }
    {
        Keyword key1 = {"BLUE", Parse_Color, 1, "main"};
        Keyword key2 = {"BLUE", Parse_Color, 1, "second"};
        if (key1==key2) ut.failure("comparison of dislike module FAILED");
    }

    {
        Keyword key = {0, Parse_Color, 1, "main"};
        if (Is_Well_Formed_Keyword(key))
            ut.failure("null keyword detect FAILED");
    }
    {
        Keyword key = {"BLUE", 0, 1, "main"};
        if (Is_Well_Formed_Keyword(key))
            ut.failure("null func detect FAILED");
    }
    {
        Keyword key = {".BLUE", Parse_Color, 1, "main"};
        if (Is_Well_Formed_Keyword(key))
            ut.failure("bad moniker detect FAILED");
    }
    {
        Keyword key = {"_BLUE", Parse_Color, 1, "main"};
        if (!Is_Well_Formed_Keyword(key))
            ut.failure("moniker with leading underscore FAILED");
    }
    {
        Keyword key = {"BLUE.", Parse_Color, 1, "main"};
        if (Is_Well_Formed_Keyword(key))
            ut.failure("bad moniker detect FAILED");
    }
}

//---------------------------------------------------------------------------//
void tstParse_Table(UnitTest &ut)
{
    Parse_Table table;
    
    table.reserve(raw_table_size);
    table.add(raw_table, raw_table_size);

    if (table.size()!=raw_table_size) ut.failure("test FAILS");
    
    File_Token_Stream token_stream("parser_test.inp");
    
    table.parse(token_stream);
    
    if (!color_set[1]) ut.failure("test FAILS");
    
    if (token_stream.error_count()!=5) ut.failure("test FAILS");


    token_stream.rewind();

    table.set_flags(Parse_Table::CASE_INSENSITIVE);

    color_set[0] = color_set[1] = 0;
    table.parse(token_stream);

    if (!color_set[1]) ut.failure("test FAILS");
    if (token_stream.error_count()!=4) ut.failure("test FAILS");

    {
        String_Token_Stream tokens("BLUE green");
        table.parse(tokens);
        if (tokens.error_count()!=0)
            ut.failure("Did NOT match mismatched case");
    }
    {
        String_Token_Stream tokens("lower blue");
        table.parse(tokens);
        if (tokens.error_count()!=0)
            ut.failure("Did NOT match lower case");
    }
    {
        String_Token_Stream tokens("lowe");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect partial match case");
    }
    {
        String_Token_Stream tokens("lower bluer");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect partial match case");
    }


    token_stream.rewind();

    table.set_flags(Parse_Table::CASE_INSENSITIVE |
		    Parse_Table::PARTIAL_IDENTIFIER_MATCH);

    color_set[0] = color_set[1] = 0;
    table.parse(token_stream);

    if (!color_set[1]) ut.failure("test FAILS");
    if (token_stream.error_count()!=3) ut.failure("test FAILS");

    // Test the Get_Flags() function, even if this test is built with the flag
    // --with-dbc=0.
    if( table.get_flags() != 3 ) ut.failure("test FAILS");
    
    // Test the check_class_invariants() function, even if this test is built
    // with the flag --with-dbc=0.
    if( ! table.check_class_invariants() ) ut.failure("test FAILS");

    // Check variations on partial match
    {
        String_Token_Stream tokens("BLUEE");
        table.parse(tokens);
        if (tokens.error_count()!=0)
            ut.failure("Did NOT match partial keyword");
    }
    {
        String_Token_Stream tokens("blue");
        table.parse(tokens);
        if (tokens.error_count()!=0)
            ut.failure("Did NOT match keyword with wrong case");
    }
    {
        String_Token_Stream tokens("end");
        if (table.parse(tokens).type()!=END)
        {
            ut.failure("END detection FAILED");
        }
    }
    // Test recovery
    {
    }
    
    table.set_flags(Parse_Table::PARTIAL_IDENTIFIER_MATCH);
    {
        String_Token_Stream tokens("BLUEE");
        table.parse(tokens);
        if (tokens.error_count()!=0)
            ut.failure("Did NOT match partial keyword");
    }
    {
        String_Token_Stream tokens("BLU green");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect mismatched case");
    }
    {
        String_Token_Stream tokens("blue");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect mismatched case");
    }
    {
        String_Token_Stream tokens("blue green red");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect mismatched case");
    }
    {
        String_Token_Stream tokens("BLUE RED");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT detect unknown keyword");
    }
    {
        String_Token_Stream tokens("BLUISH");
        table.parse(tokens);
        if (tokens.error_count()!=1)
            ut.failure("Did NOT catch partial mismatch");
    }
    {
        String_Token_Stream tokens("end");
        if (table.parse(tokens).type()!=END)
        {
            ut.failure("END detection FAILED");
        }
        if (table.parse(tokens).type()!=EXIT)
        {
            ut.failure("exit detection FAILED");
        }
    }
    {
        Colon_Token_Stream tokens;
        if (table.parse(tokens).type()!=END)
        {
            ut.failure("END detection FAILED");
        }
        if (table.parse(tokens).type()!=EXIT)
        {
            ut.failure("exit detection FAILED");
        }
    }
    // Error handling
    {
        Error_Token_Stream tokens;
        if (table.parse(tokens).type()!=ERROR)
        {
            ut.failure("error detection FAILED");
        }
    }

    Parse_Table table_2(raw_table, raw_table_size);

    if (table_2.size()!=raw_table_size) ut.failure("test FAILS");
    
    token_stream.rewind();
    
    table_2.parse(token_stream);
    
    if (!color_set[1]) ut.failure("test FAILS");
    
    if (token_stream.error_count()!=5) ut.failure("error count FAILS");

    Keyword test_key = {"THIS SHOULD WORK", Parse_Color, 0, 0};
    if (!Is_Well_Formed_Keyword(test_key)) ut.failure("test FAILS");

    Keyword benign_ambiguous_table[] = 
	{
	    {"KEY", Parse_Color, 0, 0},
	    {"KEY", Parse_Color, 0, 0}
	};
    table_2.add(benign_ambiguous_table, 2);
    token_stream.rewind();
    table_2.parse(token_stream);

    Keyword malign_ambiguous_table[] = 
	{
	    {"KEY", Parse_Color, 1, 0}
	};
    try 
    {
	table_2.add(malign_ambiguous_table, 1);
	token_stream.rewind();
	table_2.parse(token_stream);
	ut.failure("did NOT catch ambiguous keyword");
    }
    catch (invalid_argument const &msg)
    {
	cout << msg.what() << endl;
	ut.passes("successfully detected ambiguous keyword");
    }

    File_Token_Stream recover_stream("recovery.inp");
    table.parse(recover_stream);
    if (recover_stream.error_count() != 2) ut.failure("test FAILS");

    Parse_Table table_3;
    Keyword case_ambiguous_table[] = 
	{
	    {"key", Parse_Color, 0, 0},
	    {"Key", Parse_Color, 1, 0}
	};    
    try 
    {
	table_3.add(case_ambiguous_table, 2);
	table_3.parse(token_stream);
	table_3.set_flags(Parse_Table::CASE_INSENSITIVE);
	token_stream.rewind();
	table_3.parse(token_stream);
	ut.failure("did NOT catch case-dependent ambiguous keyword");
    }
    catch (invalid_argument const &msg)
    {
	cout << msg.what() << endl;
	ut.passes("successfully detected case-dependent ambiguous keyword");
    }

    Parse_Table table_3a;
    Keyword casea_ambiguous_table[] = 
	{
	    {"key", Parse_Color, 0, 0},
	    {"Key", Parse_Any_Color, 0, 0}
	};    
    try 
    {
	table_3a.add(casea_ambiguous_table, 2);
	table_3a.parse(token_stream);
	table_3a.set_flags(Parse_Table::CASE_INSENSITIVE);
	token_stream.rewind();
	table_3a.parse(token_stream);
	ut.failure("did NOT catch case-dependent ambiguous keyword");
    }
    catch (invalid_argument const &msg)
    {
	cout << msg.what() << endl;
	ut.passes("successfully detected case-dependent ambiguous keyword");
    }

    Parse_Table table_4;
    table_4.add(table_4);
    table_4.add(table);
        
    recover_stream.rewind();
    table_4.parse(recover_stream);
    if (recover_stream.error_count() != 2) ut.failure("test FAILS");

    table.set_flags(Parse_Table::ONCE);
    {
        String_Token_Stream tokens("BLUE, ERROR");
        if (table.parse(tokens).type()!=END ||
            tokens.error_count()>0)
        {
            ut.failure("FAILED to end on one token in ONCE mode");
        }
    }

    {
	// Additional test mandated by bug discovery.

	Parse_Table table;
	
	table.reserve(raw_table_2_size);
	table.add(raw_table_2, raw_table_2_size);
	
	if (table.size()!=raw_table_2_size) ut.failure("test FAILS");
	
	File_Token_Stream token_stream("parser_test.inp");
	
	table.parse(token_stream);

	if (token_stream.error_count() != 5) ut.failure("test FAILS");
    }

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
        tstKeyword(ut);
        tstParse_Table(ut);
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }

    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstParse_Table.cc
//---------------------------------------------------------------------------//
