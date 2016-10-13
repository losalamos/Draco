/*-----------------------------------*-C-*-----------------------------------*/
/*!
 * \file   parser/Debug_Options.cc
 * \author Kent Grimmett Budge
 * \brief
 * \note   Copyright (C) 2014-2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
/*---------------------------------------------------------------------------*/
/* $Id: template.h 7388 2014-01-22 16:02:07Z kellyt $ */
/*---------------------------------------------------------------------------*/

#include "Debug_Options.hh"

namespace rtt_parser {
using std::string;

//---------------------------------------------------------------------------------------//
unsigned parse_debug_options(Token_Stream &tokens) {
  using namespace rtt_parser;

  unsigned Result = 0;

  Token token = tokens.lookahead();
  while (token.type() == KEYWORD) {
    if (token.text() == "ALGORITHM") {
      Result += DEBUG_ALGORITHM;
      tokens.shift();
    } else if (token.text() == "TIMESTEP") {
      Result += DEBUG_TIMESTEP;
      tokens.shift();
    } else if (token.text() == "TIMING") {
      Result += DEBUG_TIMING;
      tokens.shift();
    } else if (token.text() == "BALANCE") {
      Result += DEBUG_BALANCE;
      tokens.shift();
    } else if (token.text() == "GMV_DUMP") {
      Result += DEBUG_GMV_DUMP;
      tokens.shift();
    } else if (token.text() == "MEMORY") {
      Result += DEBUG_MEMORY;
      tokens.shift();
    } else if (token.text() == "RESET_TIMING") {
      Result += DEBUG_RESET_TIMING;
      tokens.shift();
    } else {
      return Result;
    }
    token = tokens.lookahead();
  }
  return Result;
}

//---------------------------------------------------------------------------------------//
string debug_options_as_text(unsigned const debug_options) {
  string Result;

  if (debug_options & DEBUG_ALGORITHM) {
    Result += ", ALGORITHM";
  }
  if (debug_options & DEBUG_TIMESTEP) {
    Result += ", TIMESTEP";
  }
  if (debug_options & DEBUG_TIMING) {
    Result += ", TIMING";
  }
  if (debug_options & DEBUG_BALANCE) {
    Result += ", BALANCE";
  }
  if (debug_options & DEBUG_GMV_DUMP) {
    Result += ", GMV_DUMP";
  }
  if (debug_options & DEBUG_MEMORY) {
    Result += ", MEMORY";
  }
  if (debug_options & DEBUG_RESET_TIMING) {
    Result += ", RESET_TIMING";
  }

  return Result;
}

} // end namespace rtt_parser

/*---------------------------------------------------------------------------*/
/* end of parser/Debug_Options.cc */
/*---------------------------------------------------------------------------*/
