//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstDebug_Options.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
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
  for (unsigned i = 0; i < 33; i++) {
    string out = debug_options_as_text(i);
    String_Token_Stream tokens(out);
    unsigned j = parse_debug_options(tokens);
    ut.check(i == j, "write/read check");
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
