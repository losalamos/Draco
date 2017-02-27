//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstString_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  Unit tests for String_Token_Stream class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/String_Token_Stream.hh"
#include "parser/utilities.hh"
#include <sstream>

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstString_Token_Stream(UnitTest &ut) {
  // Build path for the input file "scanner_test.inp"
  string const stInputFile(ut.getTestSourcePath() +
                           std::string("scanner_test.inp"));

  ifstream infile(stInputFile.c_str());
  string contents;
  while (true) {
    char const c = infile.get();
    if (infile.eof() || infile.fail())
      break;
    contents += c;
  }

  {
    String_Token_Stream tokens(contents);
    if (tokens.whitespace() != Text_Token_Stream::default_whitespace)
      FAILMSG("whitespace characters are NOT correct defaults");
    else
      PASSMSG("whitespace characters are correct defaults");

    Token token = tokens.lookahead(4);
    if (token.type() != KEYWORD || token.text() != "BLACK")
      FAILMSG("lookahead(4) does NOT have correct value");
    else
      PASSMSG("lookahead(4) has correct value");

    tokens.report_semantic_error(token, "dummy error");
    if (tokens.error_count() != 1)
      FAILMSG("Dummy error NOT counted properly");
    else
      PASSMSG("Dummy error counted properly");

    if (!tokens.check_class_invariants())
      ITFAILS;
  }

  {
    set<char> ws;
    ws.insert(':');
    String_Token_Stream tokens(contents, ws);
    if (tokens.whitespace() != ws)
      FAILMSG("whitespace characters are NOT correctly specified");
    else
      PASSMSG("whitespace characters are correctly specified");

    Token token = tokens.lookahead(4);
    if (token.type() != OTHER || token.text() != "=")
      FAILMSG("lookahead(4) does NOT have correct value");
    else
      PASSMSG("lookahead(4) has correct value");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLUE")
      FAILMSG("First shift does NOT have correct value");
    else
      PASSMSG("First shift has correct value");

    token = tokens.lookahead();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      FAILMSG("Lookahed after first shift does NOT have correct value");
    else
      PASSMSG("Lookahead after first shift has correct value");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ERROR")
      FAILMSG("Second shift does NOT have correct value");
    else
      PASSMSG("Second shift has correct value");

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "GENERATE ANOTHER ERROR")
      FAILMSG("Third shift does NOT have correct value");
    else
      PASSMSG("Third shift has correct value");

    token = Token('$', "test_parser");
    tokens.pushback(token);

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "$")
      FAILMSG("Shift after pushback does NOT have correct value");
    else
      PASSMSG("Shift after pushback has correct value");

    try {
      tokens.report_syntax_error(token, "dummy syntax error");
      FAILMSG("Syntax error NOT correctly thrown");
    } catch (const Syntax_Error &msg) {
      PASSMSG("Syntax error correctly thrown and caught");
    }
    if (tokens.error_count() != 1) {
      FAILMSG("Syntax error NOT correctly counted");
    } else {
      PASSMSG("Syntax error correctly counted");
      if (tokens.messages() == "test_parser\ndummy syntax error\n")
        PASSMSG("Correct error message");
      else
        FAILMSG("NOT correct error message");
    }

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "COLOR")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "=")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLACK")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != END)
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "-")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != ".563e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != ".")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "-")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.563")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1.e3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1e+3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != REAL || token.text() != "1e3")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "19090")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "01723")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0x1111a")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "8123")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != STRING || token.text() != "\"manifest string\"")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != STRING || token.text() != "\"manifest \\\"string\\\"\"")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != OTHER || token.text() != "@")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "1")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "e")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "x")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;
    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;

    tokens.rewind();
    token = tokens.lookahead();
    token = tokens.shift();
    if (token.type() != KEYWORD || token.text() != "BLUE")
      ITFAILS;
  }

  //---------------------------------------------------------------------------//

  {

    // Build path for the input file "scanner_recovery.inp"
    string const srInputFile(ut.getTestSourcePath() +
                             std::string("scanner_recovery.inp"));

    ifstream infile(srInputFile.c_str());
    string contents;
    while (true) {
      char const c = infile.get();
      if (infile.eof() || infile.fail())
        break;
      contents += c;
    }
    String_Token_Stream tokens(contents);
    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when\n"
          << "\tunbalanced quotes were read from the input\n"
          << "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      // cout << msg.what() << endl;
      // exception = true;
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream msg;
        msg << "Caught expected exception from Token_Stream.\n"
            << "\tunbalanced quotes were read from the input\n"
            << "\tfile, \"scanner_recover.inp\" (line 1)." << endl;
        PASSMSG(msg.str());
      } else
        ITFAILS;
    }

    try {
      tokens.shift();
      ostringstream msg;
      msg << "Token_Stream did not throw an exception when\n"
          << "\tunbalanced quotes were read from the input\n"
          << "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
      FAILMSG(msg.str());
    } catch (const Syntax_Error &msg) {
      //cout << msg.what() << endl;
      // exception = true;
      string errmsg = msg.what();
      string expected("syntax error");
      if (errmsg == expected) {
        ostringstream msg;
        msg << "Caught expected exception from Token_Stream.\n"
            << "\tunbalanced quotes were read from the input\n"
            << "\tfile, \"scanner_recover.inp\" (line 2)." << endl;
        PASSMSG(msg.str());
      } else
        ITFAILS;
    }

    // Now test assignment of new string and dipthong OTHER tokens

    tokens = string("<= >= && ||");

    Token token = tokens.shift();
    if (token.text() != "<=")
      ITFAILS;

    token = tokens.shift();
    if (token.text() != ">=")
      ITFAILS;

    token = tokens.shift();
    if (token.text() != "&&")
      ITFAILS;

    token = tokens.shift();
    if (token.text() != "||")
      ITFAILS;

    token = tokens.shift();
    if (token.type() != EXIT)
      ITFAILS;
  }

  {
    String_Token_Stream tokens("09");
    if (tokens.is_nb_whitespace('\t')) {
      ut.passes("tab correctly identified as nonbreaking whitespace");
    } else {
      ut.failure("tab NOT correctly identified as nonbreaking "
                 "whitespace");
    }
    Token token = tokens.shift();
    if (token.type() != INTEGER || token.text() != "0") {
      ut.failure("did NOT scan 09 correctly");
    } else {
      ut.passes("scanned 09 correctly");
    }
  }
  {
    String_Token_Stream tokens("_, __, _ _, > < & | 1E3 0XA");
    if (tokens.shift().text() != "_")
      ut.failure("Did NOT correctly scan _");
    if (tokens.shift().text() != "__")
      ut.failure("Did NOT correctly scan __");
    if (tokens.shift().text() != "_ _")
      ut.failure("Did NOT correctly scan _ _");
    if (tokens.shift().text() != ">")
      ut.failure("Did NOT correctly scan >");
    if (tokens.shift().text() != "<")
      ut.failure("Did NOT correctly scan <");
    if (tokens.shift().text() != "&")
      ut.failure("Did NOT correctly scan &");
    if (tokens.shift().text() != "|")
      ut.failure("Did NOT correctly scan |");
    if (!soft_equiv(parse_real(tokens), 1e3))
      ut.failure("Did NOT correctly scan 1E3");
    if (parse_integer(tokens) != 10)
      ut.failure("Did NOT correctly scan 0XA");
  }

  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstString_Token_Stream(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstString_Token_Stream.cc
//---------------------------------------------------------------------------//
