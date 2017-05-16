//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstClass_Parser.cc
 * \author Kent Budge
 * \date   Mon Aug 28 07:36:50 2006
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Class_Parse_Table.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/utilities.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
class DummyClass {
public:
  DummyClass(double const insouciance) : insouciance(insouciance) {}

  double Get_Insouciance() const { return insouciance; }

private:
  double insouciance;
};

namespace rtt_parser {
//---------------------------------------------------------------------------//
template <> class Class_Parse_Table<DummyClass> {
public:
  // TYPEDEFS

  typedef DummyClass Return_Class;

  // MANAGEMENT

  Class_Parse_Table(bool is_indolent = false);

  // SERVICES

  Parse_Table const &parse_table() const { return parse_table_; }

  bool allow_exit() const { return false; }

  void check_completeness(Token_Stream &tokens);

  std::shared_ptr<DummyClass> create_object();

protected:
  // DATA

  double parsed_insouciance;

  bool is_indolent;

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
template <>
std::shared_ptr<DummyClass> parse_class<DummyClass>(Token_Stream &tokens) {
  return parse_class_from_table<Class_Parse_Table<DummyClass>>(tokens);
}

//---------------------------------------------------------------------------//
template <>
std::shared_ptr<DummyClass> parse_class<DummyClass>(Token_Stream &tokens,
                                                    bool const &is_indolent) {
  return parse_class_from_table<Class_Parse_Table<DummyClass>>(tokens,
                                                               is_indolent);
}

//---------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::parse_insouciance_(Token_Stream &tokens,
                                                       int) {
  tokens.check_semantics(current_->parsed_insouciance < 0.0,
                         "duplicate specification of insouciance");

  current_->parsed_insouciance = parse_real(tokens);
  if (current_->parsed_insouciance < 0) {
    tokens.report_semantic_error("insouciance must be nonnegative");
    current_->parsed_insouciance = 1;
  }

  if (current_->is_indolent) {
    current_->parsed_insouciance = -current_->parsed_insouciance;
  }
}

//---------------------------------------------------------------------------//
Class_Parse_Table<DummyClass>::Class_Parse_Table(bool const is_indolent)
    : parsed_insouciance(-1.0) // sentinel value
{
  if (!parse_table_is_initialized_) {
    Keyword const keywords[] = {
        {"insouciance", parse_insouciance_, 0, ""},
    };
    unsigned const number_of_keywords = sizeof(keywords) / sizeof(Keyword);

    parse_table_.add(keywords, number_of_keywords);

    parse_table_is_initialized_ = true;
  }

  this->is_indolent = is_indolent;

  current_ = this;
}

//---------------------------------------------------------------------------//
void Class_Parse_Table<DummyClass>::check_completeness(Token_Stream &tokens) {
  tokens.check_semantics(is_indolent || parsed_insouciance >= 0,
                         "insouciance was not specified");
}

//---------------------------------------------------------------------------//
std::shared_ptr<DummyClass> Class_Parse_Table<DummyClass>::create_object() {
  std::shared_ptr<DummyClass> Result =
      std::shared_ptr<DummyClass>(new DummyClass(parsed_insouciance));
  return Result;
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstClass_Parser(UnitTest &ut) {
  string text = "insouciance = 3.3\nend\n";
  String_Token_Stream tokens(text);

  std::shared_ptr<DummyClass> dummy = parse_class<DummyClass>(tokens);

  ut.check(dummy != nullptr, "parsed the class object", true);
  ut.check(dummy->Get_Insouciance() == 3.3, "parsed the insouciance correctly");

  tokens.rewind();
  dummy = parse_class<DummyClass>(tokens, true);

  ut.check(dummy != nullptr, "parsed the indolent class object", true);
  ut.check(dummy->Get_Insouciance() == -3.3,
           "parsed the indolent insouciance correctly");

  // Test that missing end is caught.

  text = "insouciance = 3.3\n";
  String_Token_Stream etokens(text);

  bool good = false;
  try {
    std::shared_ptr<DummyClass> dummy = parse_class<DummyClass>(etokens);
  } catch (Syntax_Error &) {
    good = true;
  }
  ut.check(good, "catches missing end");

  tokens.rewind();
  good = false;
  try {
    std::shared_ptr<DummyClass> dummy = parse_class<DummyClass>(etokens, true);
  } catch (Syntax_Error &) {
    good = true;
  }
  ut.check(good, "indolent catches missing end");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstClass_Parser(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstClass_Parser.cc
//---------------------------------------------------------------------------//
