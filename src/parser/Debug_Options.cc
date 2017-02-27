/*-----------------------------------*-C-*-----------------------------------*/
/*!
 * \file   parser/Debug_Options.cc
 * \author Kent Grimmett Budge
 * \brief
 * \note   Copyright (C) 2014-2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
/*---------------------------------------------------------------------------*/
/* $Id: template.h 7388 2014-01-22 16:02:07Z kellyt $ */
/*---------------------------------------------------------------------------*/

#include "Debug_Options.hh"

namespace rtt_parser {
using std::string;

//---------------------------------------------------------------------------------------//
/*! Parse a debug specification.
 *
 * \param tokens Token stream from which to parse a debug specification. The specification
 * is a set of debug keywords, each optionally prefixed with a '!', and ends with the first
 * token that is not a recognized debug keyword.
 *
 * \param parent Optional parent mask; defaults to zero. Allows a debug mask to be
 * based on a parent mask, with selected bits added or masked out.
 *
 * \return A debug options mask.
 */
unsigned parse_debug_options(Token_Stream &tokens, unsigned const parent) {
  using namespace rtt_parser;

  unsigned Result = parent;

  Token token = tokens.lookahead();
  while (token.type() == KEYWORD || token.text() == "!") {
    bool mask_out = (token.text() == "!");
    if (mask_out) {
      tokens.shift();
      token = tokens.lookahead();
    }
    unsigned mask = 0;
    if (token.text() == "ALGORITHM") {
      mask = DEBUG_ALGORITHM;
    } else if (token.text() == "TIMESTEP") {
      mask = DEBUG_TIMESTEP;
    } else if (token.text() == "TIMING") {
      mask = DEBUG_TIMING;
    } else if (token.text() == "BALANCE") {
      mask = DEBUG_BALANCE;
    } else if (token.text() == "GMV_DUMP") {
      mask = DEBUG_GMV_DUMP;
    } else if (token.text() == "MEMORY") {
      mask = DEBUG_MEMORY;
    } else if (token.text() == "RESET_TIMING") {
      mask = DEBUG_RESET_TIMING;
    }
    if (mask) {
      if (mask_out) {
        Result = Result & ~mask;
      } else {
        Result = Result | mask;
      }
      tokens.shift();
    } else {
      tokens.check_syntax(!mask_out, "trailing '!'");
      return Result;
    }
    token = tokens.lookahead();
  }
  return Result;
}

//---------------------------------------------------------------------------------------//
/*! Convert a debug mask to a string containing comma-delimited set of debug keywords.
 *
 * \param debug_options Debug mask to be converted to a set of keywords.
 *
 * \return A string containing a comma-delimited set of debug options.
 */
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
