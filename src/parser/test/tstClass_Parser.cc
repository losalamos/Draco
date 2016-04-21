//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstClass_Parser.cc
 * \author Kent Budge
 * \date   Mon Aug 28 07:36:50 2006
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/utilities.hh"
#include "parser/Class_Parse_Table.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
class DummyClass
{
  public:

    DummyClass(double const insouciance)
        :
        insouciance(insouciance)
    {
        Require(insouciance>=0.0);
    }

    double Get_Insouciance() const { return insouciance; }

  private:

    double insouciance;
};

namespace rtt_parser
{
//---------------------------------------------------------------------------//
template<>
class Class_Parse_Table<DummyClass>
{
  public:

    // TYPEDEFS

    typedef DummyClass Return_Class;

    // MANAGEMENT

    Class_Parse_Table();

    // SERVICES

    Parse_Table const &parse_table() const { return parse_table_; }

    bool allow_exit() const { return false; }

    void check_completeness(Token_Stream &tokens);

    SP<DummyClass> create_object();

  protected:

    // DATA

    double parsed_insouciance;

  private:

    // IMPLEMENTATION

    static void parse_insouciance_(Token_Stream &tokens, int);

    // STATIC

    static Class_Parse_Table *current_;
    static Parse_Table parse_table_;
    static bool parse_table_is_initialized_;
};

//---------------------------------------------------------------------------//
Class_Parse_Table<DummyClass> *Class_Parse_Table<DummyClass>::current_;
Parse_Table Class_Parse_Table<DummyClass>::parse_table_;
bool Class_Parse_Table<DummyClass>::parse_table_is_initialized_ = false;

//---------------------------------------------------------------------------//
template<>
SP<DummyClass> parse_class<DummyClass>(Token_Stream &tokens)
{
    return parse_class_from_table<Class_Parse_Table<DummyClass> >(tokens);
}

//---------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::parse_insouciance_(Token_Stream &tokens, int)
{
    if (current_->parsed_insouciance>=0.0)
    {
        tokens.report_semantic_error("duplicate specification of insouciance");
    }
    current_->parsed_insouciance = parse_real(tokens);
    if (current_->parsed_insouciance<0)
    {
        tokens.report_semantic_error("insouciance must be nonnegative");
        current_->parsed_insouciance = 1;
    }
}

//---------------------------------------------------------------------------//
Class_Parse_Table<DummyClass>::Class_Parse_Table()
    :  parsed_insouciance(-1.0) // sentinel value
{
    if (!parse_table_is_initialized_)
    {
        Keyword const keywords[] =
            {
                {"insouciance", parse_insouciance_, 0, ""},
            };
        unsigned const number_of_keywords = sizeof(keywords)/sizeof(Keyword);

        parse_table_.add(keywords, number_of_keywords);

        parse_table_is_initialized_ = true;
    }
    current_ = this;
}

//---------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::check_completeness(Token_Stream &tokens)
{
    if (parsed_insouciance<0)
    {
        tokens.report_semantic_error("insouciance was not specified");
    }
}

//---------------------------------------------------------------------------//
SP<DummyClass> Class_Parse_Table<DummyClass>::create_object()
{
    SP<DummyClass> Result = SP<DummyClass>(new DummyClass(parsed_insouciance));
    return Result;
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstClass_Parser( UnitTest & ut)
{
    string text = "insouciance = 3.3\nend\n";
    String_Token_Stream tokens(text);

    SP<DummyClass> dummy = parse_class<DummyClass>(tokens);

    if (dummy != SP<DummyClass>())
    {
        ut.passes("parsed the class object");

        if (dummy->Get_Insouciance() == 3.3)
        {
            ut.passes("parsed the insouciance correctly");
        }
        else
        {
            ut.failure("did NOT parse the insouciance correctly");
        }
    }
    else
    {
        cout << tokens.messages() << endl;
        ut.failure("did NOT parse the class object");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try { tstClass_Parser(ut); }
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstClass_Parser.cc
//---------------------------------------------------------------------------//
