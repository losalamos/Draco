//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstDebug_Options.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstDebug_Options.cc 8318 2016-04-21 03:04:14Z kellyt $
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Debug_Options.hh"
#include "parser/String_Token_Stream.hh"

using namespace std;

using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void debug_options_test(UnitTest &ut) {
  for (unsigned i = 1; i < 2 * DEBUG_RESET_TIMING - 1; i++) {
    string out = debug_options_as_text(i);
    String_Token_Stream tokens(out);
    unsigned j = parse_debug_options(tokens);
    ut.check(i == j, "write/read check");
  }
  for (unsigned i = 1; i <= DEBUG_RESET_TIMING; i <<= 1U) {
    string out = debug_options_as_text(i);
    out = '!' + out.substr(1);
    String_Token_Stream mask_tokens(out);
    unsigned j = parse_debug_options(mask_tokens, i);
    ut.check(j == 0, "write/read mask check");
  }
  bool did = true;
  try {
    string out = "!";
    String_Token_Stream tokens(out);
    parse_debug_options(tokens);
  } catch (Syntax_Error &) {
    did = false;
  }
  ut.check(!did, "catches syntax error for trailing '!'");

  // try adding a new keyword then removing it
  {
    unsigned bit = add_debug_option("EXTENDED");
    // duplicate okay
    ut.check(add_debug_option("EXTENDED") == bit,
             "duplicate definition benign");
    string out = "EXTENDED";
    String_Token_Stream tokens(out);
    ut.check(parse_debug_options(tokens) == bit, "added debug option");
    flush_debug_options();
    ut.check(get_debug_option("EXTENDED") == 0, "flushed debug option");
    // be sure we didn't flush standard options
    ut.check(get_debug_option("TIMING") != 0, "did not flush standard");
  }

  // try assigning a new keyword
  {
    add_debug_option("EXTENDED", DEBUG_END);
    add_debug_option("EXTENDED", DEBUG_END); // should be benign
    string out = "EXTENDED";
    String_Token_Stream tokens(out);
    ut.check(parse_debug_options(tokens) == DEBUG_END, "added debug option");
    bool did = true;
    try {
      add_debug_option("EXTENDED", 2 * DEBUG_END); // duplicate inconsistent
    } catch (std::invalid_argument &) {
      did = false;
    }
    ut.check(!did, "catches illegal argument error for inconsistent duplicate");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    debug_options_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstDebug_Options.cc
//---------------------------------------------------------------------------//
