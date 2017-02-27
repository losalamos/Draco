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
